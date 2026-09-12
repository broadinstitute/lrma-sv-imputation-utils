#!/usr/bin/env python3
"""Generate fixtures + expected output for pop-glimpse2-joint-opt.

The tool does f32 arithmetic, so this is a *faithful f32 re-implementation*
(every intermediate is rounded to f32 via struct pack/unpack, exactly as the
Rust binary does). Because a handful of steps round to a grid --

  * GT hap call at the 0.5 probability threshold,
  * DS  = format_float(p0+p1)      -> nearest 0.001,
  * GP  = round(gp_raw * 1000)     -> nearest permille,
  * top-k allele cut when max_alleles < num_alleles,

-- we use REJECTION SAMPLING: random (seeded) GP inputs are drawn, the oracle
computes every boundary-sensitive quantity, and any draw that lands closer than
a comfortable margin to *any* boundary is rejected. The surviving fixture is
therefore guaranteed to produce byte-identical output from the Rust f32 code,
independent of tie-breaking or half-even-vs-half-away rounding differences.

Two expected files are produced from the SAME fixtures, one per max_alleles
value, so the top-k path is covered without a second fixture set.
"""

import os
import random
import struct
import sys

# --------------------------------------------------------------------------- #
# f32 helpers -- mirror Rust's `as f32` / f32 arithmetic exactly.
# --------------------------------------------------------------------------- #
def f32(x):
    return struct.unpack("f", struct.pack("f", float(x)))[0]

F32_1 = f32(1.0)
CLAMP_LO = f32(1e-5)
CLAMP_HI = f32(1.0 - 1e-5)      # computed in f32, matching `1.0 - 1e-5`


def clamp(v, lo, hi):
    if v < lo:
        return lo
    if v > hi:
        return hi
    return v


def rha(x):
    """round-half-away-from-zero for non-negative x, like Rust f32::round."""
    import math
    return math.floor(x + 0.5)


def format_float(val):
    """Port of the tool's format_float: `{:.3}` then strip trailing 0 / '.'."""
    s = "%.3f" % val
    s = s.rstrip("0").rstrip(".")
    return s if s else "0"


# --------------------------------------------------------------------------- #
# Faithful port of process_group. Returns (body_lines, min_slack) where
# min_slack is the smallest distance to any rounding / threshold boundary
# encountered (used only by the generator's rejection loop).
# --------------------------------------------------------------------------- #
def process_group(group, id_buffer, max_alleles):
    """group: list of dicts, each an allele line:
         {info: {RAF,AF,INFO}, atomic_ids: [..], samples: [(gt_str, gp_triple)]}
       id_buffer: key -> (pos, orig_id, ref, alt, ac, an)
    """
    slack = float("inf")
    num_alleles = len(group)
    num_samples = len(group[0]["samples"])
    chrom = group[0]["chrom"]
    
    max_bubble_info = -1.0
    has_max_bubble = False
    for rec in group:
        if rec["info"].get("INFO") is not None:
            max_bubble_info = max(max_bubble_info, float(rec["info"]["INFO"]))
            has_max_bubble = True

    all_atomic = []
    seen = set()
    # hap_probs[s][a] = (p0, p1)
    hap_probs = [[(0.0, 0.0)] * num_alleles for _ in range(num_samples)]

    for a, rec in enumerate(group):
        for aid in rec["atomic_ids"]:
            if aid not in seen:
                seen.add(aid)
                all_atomic.append(aid)
        for s in range(num_samples):
            gt_val, (g0, g1, g2) = rec["samples"][s]
            gp1 = f32(g1)
            gp2 = f32(g2)
            if gt_val.startswith("1|0"):
                p0, p1 = f32(gp2 + gp1), gp2
            elif gt_val.startswith("0|1"):
                p0, p1 = gp2, f32(gp2 + gp1)
            else:
                half = f32(gp2 + f32(gp1 / f32(2.0)))
                p0, p1 = half, half
            hap_probs[s][a] = (clamp(p0, CLAMP_LO, CLAMP_HI),
                               clamp(p1, CLAMP_LO, CLAMP_HI))

    # accumulate per-atomic distributions
    dists = {aid: [[0.0, 0.0] for _ in range(num_samples)] for aid in all_atomic}

    for s in range(num_samples):
        scores = []
        for a in range(num_alleles):
            scores.append((a, f32(hap_probs[s][a][0] + hap_probs[s][a][1])))
        # stable descending sort by f32 score (ties keep ascending index)
        order = sorted(range(num_alleles), key=lambda a: scores[a][1], reverse=True)
        m = min(num_alleles, max_alleles)
        if m < num_alleles:
            cut = scores[order[m - 1]][1] - scores[order[m]][1]
            slack = min(slack, abs(cut))          # score-cut margin
        top = order[:m]

        z0, z1 = F32_1, F32_1
        w0m, w1m = {}, {}
        for a in top:
            p0, p1 = hap_probs[s][a]
            w0 = f32(p0 / f32(F32_1 - p0))
            w1 = f32(p1 / f32(F32_1 - p1))
            w0m[a] = w0
            w1m[a] = w1
            z0 = f32(z0 + w0)
            z1 = f32(z1 + w1)
        for a in top:
            np0 = f32(w0m[a] / z0)
            np1 = f32(w1m[a] / z1)
            for aid in group[a]["atomic_ids"]:
                d = dists[aid][s]
                d[0] = f32(d[0] + np0)
                d[1] = f32(d[1] + np1)

    # sort atomic vars by (coord, ref, alt) from id_buffer
    sorted_vars = sorted(all_atomic,
                         key=lambda a: (id_buffer[a][0], id_buffer[a][2], id_buffer[a][3]))

    lines = []
    for aid in sorted_vars:
        pos, orig_id, ref, alt, ac, an = id_buffer[aid]
        
        info = ["ID=%s" % aid]
        if an > 0:
            info.append("RAF=%s" % format_float(f32(ac / an)))

        cols = [chrom, str(pos), orig_id, ref, alt, ".", ".", "INFO_PH", "GT:DS:GP"]
        
        ds_sum = 0.0
        ds2_sum = 0.0
        ds4_sum = 0.0
        sample_strings = []

        for s in range(num_samples):
            p0 = clamp(dists[aid][s][0], 0.0, F32_1)
            p1 = clamp(dists[aid][s][1], 0.0, F32_1)

            slack = min(slack, abs(p0 - 0.5), abs(p1 - 0.5))   # GT threshold
            hap0 = "1" if p0 > 0.5 else "0"
            hap1 = "1" if p1 > 0.5 else "0"

            ds_raw = f32(p0 + p1)
            slack = min(slack, _round_slack(ds_raw * 1000.0))  # DS 0.001 grid
            ds = format_float(ds_raw)

            gp0 = f32(f32(F32_1 - p0) * f32(F32_1 - p1))
            gp1 = f32(f32(p0 * f32(F32_1 - p1)) + f32(f32(F32_1 - p0) * p1))
            gp2 = f32(p0 * p1)
            for g in (gp0, gp1, gp2):
                slack = min(slack, _round_slack(g * 1000.0))   # permille grid

            v0 = int(rha(f32(gp0 * 1000.0)))
            v1 = int(rha(f32(gp1 * 1000.0)))
            v2 = int(rha(f32(gp2 * 1000.0)))
            diff = 1000 - (v0 + v1 + v2)
            if diff != 0:
                if v0 >= v1 and v0 >= v2:
                    v0 += diff
                elif v1 >= v0 and v1 >= v2:
                    v1 += diff
                else:
                    v2 += diff
            
            gp1_q = f32(v1 / 1000.0)
            gp2_q = f32(v2 / 1000.0)
            ds_q = f32(gp1_q + f32(2.0 * gp2_q))
            ds_sum = f32(ds_sum + ds_q)
            ds2_sum = f32(ds2_sum + f32(ds_q * ds_q))
            ds4_sum = f32(ds4_sum + f32(gp1_q + f32(4.0 * gp2_q)))
            
            gp = "%s,%s,%s" % (format_float(f32(v0 / 1000.0)),
                               format_float(gp1_q),
                               format_float(gp2_q))
            sample_strings.append("%s:%s:%s" % ("%s|%s" % (hap0, hap1), ds, gp))

        n_tar_haps = f32(2.0 * num_samples)
        safe_n = max(n_tar_haps, f32(1e-9))
        af = f32(ds_sum / safe_n)
        denom = f32(n_tar_haps * f32(af * f32(F32_1 - af)))

        recalc_info = F32_1
        if af > 0.0 and af < 1.0 and denom > 0.0:
            recalc_info = f32(F32_1 - f32(f32(ds4_sum - ds2_sum) / denom))
        recalc_info = max(recalc_info, 0.0)
        recalc_info_rounded = f32(rha(f32(recalc_info * 1000.0)) / 1000.0)

        info.append("AF=%s" % format_float(af))
        info.append("INFO=%s" % format_float(recalc_info_rounded))
        if has_max_bubble:
            info.append("INFO_MAX_BUBBLE=%s" % format_float(max_bubble_info))

        cols[7] = ";".join(info)
        cols.extend(sample_strings)
        lines.append("\t".join(cols))
    return lines, slack


def _round_slack(permille_value):
    """distance to the nearest .5 rounding tie, in permille units."""
    frac = permille_value - int(permille_value)
    if frac < 0:
        frac += 1.0
    return abs(frac - 0.5)


# --------------------------------------------------------------------------- #
# Fixture scenario. Deterministic structure; only the GP numbers are sampled.
# --------------------------------------------------------------------------- #
NUM_SAMPLES = 3
MAX_ALLELES_CASES = [10, 2]     # >=num_alleles (full) and <num_alleles (top-k)

# atomic variants: key -> (pos, orig_id, ref, alt, ac, an)
ID_BUFFER = {
    "v1": (1000, "rs1", "A", "G", 100, 1000),
    "v2": (1000, "rs2", "A", "T", 200, 1000),
    "v3": (1001, "rs3", "C", "T", 300, 1000),
    "v4": (2000, "rs4", "G", "A", 400, 1000),
    "v5": (2000, "rs5", "G", "C", 500, 1000),
}

# One bubble at chr1:1000 (3 alt lines) and one at chr1:2000 (3 alt lines).
# Each entry: (chrom, pos, ref, alt, atomic_ids, has_raf, has_af, has_info)
BUBBLES = [
    ("chr1", 1000, [
        ("A", "G", ["v1"],       True,  True,  True),
        ("A", "T", ["v1", "v2"], True,  False, True),   # shares v1 -> accumulation
        ("A", "C", ["v3"],       False, True,  False),
    ]),
    ("chr1", 2000, [
        ("G", "A", ["v4"],       True,  True,  False),
        ("G", "C", ["v5"],       True,  True,  True),
        ("G", "T", ["v4", "v5"], False, False, True),   # shares both
    ]),
]

GT_CHOICES = ["1|0", "0|1", "1|1", "0|0"]


def sample_gp(rng):
    """A clean permille GP triple that sums to 1.0 (values as .3f strings)."""
    a = rng.randint(20, 960)
    b = rng.randint(20, 1000 - a - 20)
    c = 1000 - a - b
    return ("%.3f" % (a / 1000.0), "%.3f" % (b / 1000.0), "%.3f" % (c / 1000.0))


def build_group(bubble, rng):
    chrom, pos, alts = bubble
    group = []
    for (ref, alt, aids, hr, ha, hi) in alts:
        samples = []
        for _ in range(NUM_SAMPLES):
            gt = rng.choice(GT_CHOICES)
            gp = sample_gp(rng)
            samples.append((gt, (float(gp[0]), float(gp[1]), float(gp[2]))))
        info = {}
        info["RAF"] = ("%.4f" % rng.uniform(0.01, 0.99)) if hr else None
        info["AF"] = ("%.4f" % rng.uniform(0.01, 0.99)) if ha else None
        info["INFO"] = ("%.3f" % rng.uniform(0.1, 0.99)) if hi else None
        group.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "atomic_ids": aids, "info": info, "samples": samples,
        })
    return group


MARGIN = 0.03   # prob units for GT; permille-frac units for rounding grids


def generate(seed=20240717):
    rng = random.Random(seed)
    while True:
        groups = [build_group(b, rng) for b in BUBBLES]
        ok = True
        for max_alleles in MAX_ALLELES_CASES:
            for g in groups:
                _, slack = process_group(g, ID_BUFFER, max_alleles)
                if slack < MARGIN:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return groups


# --------------------------------------------------------------------------- #
# File writers
# --------------------------------------------------------------------------- #
CONTIG = "##contig=<ID=chr1,length=100000>"
SAMPLE_NAMES = ["S%d" % (i + 1) for i in range(NUM_SAMPLES)]

# Header lines placed in the MAIN vcf: the AK/GL/KC ones must be dropped by the
# tool, the rest passed through verbatim.
MAIN_HEADER = [
    "##fileformat=VCFv4.2",
    CONTIG,
    '##INFO=<ID=AK,Number=1,Type=String,Description="dropped by tool">',
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="dropped by tool">',
    '##FORMAT=<ID=KC,Number=1,Type=Integer,Description="dropped by tool">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">',
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="genotype probs">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLE_NAMES),
]
DROP_TOKENS = ("INFO=<ID=AK", "FORMAT=<ID=GL", "FORMAT=<ID=KC")


def info_str(info):
    parts = []
    if info.get("RAF") is not None:
        parts.append("RAF=%s" % info["RAF"])
    if info.get("AF") is not None:
        parts.append("AF=%s" % info["AF"])
    if info.get("INFO") is not None:
        parts.append("INFO=%s" % info["INFO"])
    return ";".join(parts) if parts else "."


def write_main(path, groups):
    with open(path, "w") as f:
        for h in MAIN_HEADER:
            f.write(h + "\n")
        for g in groups:
            for rec in g:
                cells = []
                for (gt, gp) in rec["samples"]:
                    cells.append("%s:%s,%s,%s" % (gt, "%.3f" % gp[0],
                                                  "%.3f" % gp[1], "%.3f" % gp[2]))
                f.write("\t".join([rec["chrom"], str(rec["pos"]), ".",
                                   rec["ref"], rec["alt"], ".", ".",
                                   info_str(rec["info"]), "GT:GP"] + cells) + "\n")


def write_sites(path, groups):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(CONTIG + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for g in groups:
            for rec in g:
                ids = rec["atomic_ids"]
                # exercise both separators the tool accepts (',' -> ':')
                sep = "," if len(ids) <= 1 else ":"
                idfield = "ID=" + sep.join(ids)
                f.write("\t".join([rec["chrom"], str(rec["pos"]), ".",
                                   rec["ref"], rec["alt"], ".", ".", idfield]) + "\n")


def write_ids(path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(CONTIG + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key, (pos, orig, ref, alt, ac, an) in sorted(ID_BUFFER.items(),
                                                         key=lambda kv: (kv[1][0], kv[0])):
            f.write("\t".join(["chr1", str(pos), orig, ref, alt, ".", ".",
                               "ID=%s;AC=%d;AN=%d" % (key, ac, an)]) + "\n")


def expected_body(groups, max_alleles):
    out = []
    for g in groups:
        lines, _ = process_group(g, ID_BUFFER, max_alleles)
        out.extend(lines)
    return out


def write_expected(path, groups, max_alleles):
    with open(path, "w") as f:
        # tool passes through main header lines except the dropped tokens
        for h in MAIN_HEADER:
            if not any(tok in h for tok in DROP_TOKENS):
                if h.startswith("#CHROM"):
                    f.write('##INFO=<ID=RAF,Number=A,Type=Float,Description="Panel reference allele frequency">\n')
                    f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Recalculated allele frequency">\n')
                    f.write('##INFO=<ID=INFO,Number=A,Type=Float,Description="Recalculated IMPUTE INFO score">\n')
                    f.write('##INFO=<ID=INFO_MAX_BUBBLE,Number=A,Type=Float,Description="Maximum INFO score across the parent bubble">\n')
                f.write(h + "\n")
        for line in expected_body(groups, max_alleles):
            f.write(line + "\n")


def main():
    here = os.environ.get("GEN_OUT_DIR") or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(here, exist_ok=True)
    exp = os.path.join(here, "expected")
    os.makedirs(exp, exist_ok=True)

    groups = generate()

    write_main(os.path.join(here, "main.vcf"), groups)
    write_sites(os.path.join(here, "sites.vcf"), groups)
    write_ids(os.path.join(here, "ids.vcf"))

    write_expected(os.path.join(exp, "max10.txt"), groups, 10)
    write_expected(os.path.join(exp, "max2.txt"), groups, 2)

    print("pop-glimpse2 fixtures + expected written "
          "(seed-stable, margin>=%.3f from all boundaries)" % MARGIN)


if __name__ == "__main__":
    main()
