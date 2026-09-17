#!/usr/bin/env python3
"""Generate fixtures + expected output for pop-glimpse2.

The tool does f32 arithmetic, so this is a *faithful f32 re-implementation*
(every intermediate is rounded to f32 via struct pack/unpack, exactly as the
Rust binary does).
"""

import os
import random
import struct
import math

def f32(x):
    return struct.unpack("f", struct.pack("f", float(x)))[0]

F32_1 = f32(1.0)
CLAMP_LO = f32(1e-5)
CLAMP_HI = f32(1.0 - 1e-5)

def clamp(v, lo, hi):
    if v < lo: return lo
    if v > hi: return hi
    return v

def rha(x):
    return math.floor(x + 0.5)

def format_float(val):
    s = "%.3f" % val
    s = s.rstrip("0").rstrip(".")
    return s if s else "0"

def process_group(group, id_buffer, max_alleles):
    slack = float("inf")
    num_alleles = len(group)
    num_samples = len(group[0]["samples"])
    chrom = group[0]["chrom"]

    all_atomic = []
    seen = set()
    carriers = {}
    hap_probs = [[(0.0, 0.0)] * num_alleles for _ in range(num_samples)]

    for a, rec in enumerate(group):
        for aid in rec["atomic_ids"]:
            if aid not in seen:
                seen.add(aid)
                all_atomic.append(aid)
            carriers.setdefault(aid, []).append(a)
        
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

    dists = {aid: [[0.0, 0.0] for _ in range(num_samples)] for aid in all_atomic}

    for s in range(num_samples):
        scores = []
        for a in range(num_alleles):
            scores.append((a, f32(hap_probs[s][a][0] + hap_probs[s][a][1])))
        
        order = sorted(range(num_alleles), key=lambda a: scores[a][1], reverse=True)
        m = min(num_alleles, max_alleles)
        if m < num_alleles:
            cut = scores[order[m - 1]][1] - scores[order[m]][1]
            slack = min(slack, abs(cut))          
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

    sorted_vars = sorted(all_atomic, key=lambda a: (id_buffer[a][0], id_buffer[a][2], id_buffer[a][3]))

    lines = []
    for aid in sorted_vars:
        pos, orig_id, ref, alt, raf = id_buffer[aid]
        carrying = carriers[aid]
        
        s_ds, s_var = 0.0, 0.0
        for s in range(num_samples):
            q0 = clamp(dists[aid][s][0], 0.0, 1.0)
            if q0 < 2e-4: q0 = 0.0
            q1 = clamp(dists[aid][s][1], 0.0, 1.0)
            if q1 < 2e-4: q1 = 0.0
            s_ds += q0 + q1
            s_var += q0 * (1.0 - q0) + q1 * (1.0 - q1)

        n_haps = 2.0 * num_samples
        af = s_ds / n_haps
        
        if 0.0 < af < 1.0:
            info_val = clamp(1.0 - s_var / (n_haps * af * (1.0 - af)), 0.0, 1.0)
        else:
            info_val = 1.0

        info = ["ID=%s" % aid]
        if raf:
            info.append("RAF=%s" % raf)
        info.append("AF=%.6f" % af)
        info.append("INFO=%.3f" % info_val)
        info.append(f"N_PATHS={len(carrying)};N_PATHS_TOTAL={num_alleles}")

        cols = [chrom, str(pos), orig_id, ref, alt, ".", ".", ";".join(info), "GT:DS:GP"]

        for s in range(num_samples):
            p0 = clamp(dists[aid][s][0], 0.0, F32_1)
            p1 = clamp(dists[aid][s][1], 0.0, F32_1)

            slack = min(slack, abs(p0 - 0.5), abs(p1 - 0.5))   
            hap0 = "1" if p0 > 0.5 else "0"
            hap1 = "1" if p1 > 0.5 else "0"

            ds_raw = f32(p0 + p1)
            slack = min(slack, _round_slack(ds_raw * 1000.0))  
            ds = format_float(ds_raw)

            gp0 = f32(f32(F32_1 - p0) * f32(F32_1 - p1))
            gp1 = f32(f32(p0 * f32(F32_1 - p1)) + f32(f32(F32_1 - p0) * p1))
            gp2 = f32(p0 * p1)
            for g in (gp0, gp1, gp2):
                slack = min(slack, _round_slack(g * 1000.0))   

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
            gp = "%s,%s,%s" % (format_float(f32(v0 / 1000.0)),
                               format_float(f32(v1 / 1000.0)),
                               format_float(f32(v2 / 1000.0)))
            cols.append("%s:%s:%s" % ("%s|%s" % (hap0, hap1), ds, gp))
        lines.append("\t".join(cols))
    return lines, slack

def _round_slack(permille_value):
    frac = permille_value - int(permille_value)
    if frac < 0:
        frac += 1.0
    return abs(frac - 0.5)

NUM_SAMPLES = 3
MAX_ALLELES_CASES = [10, 2, 1]     

ID_BUFFER = {
    "v1": (1000, "rs1", "A", "G", "0.200000"),
    "v2": (1000, "rs2", "A", "T", "0.150000"),
    "v3": (1001, "rs3", "C", "T", "0.050000"),
    "v4": (2000, "rs4", "G", "A", "0.300000"),
    "v5": (2000, "rs5", "G", "C", "0.400000"),
    "v6": (3000, "rs6", "T", "C", "0.000000"),
    "v7": (4000, "rs7", "C", "A", "0.010000"),
    "v8": (5000, "rs8", "A", "G", "0.500000"),
    "v9": (5001, "rs9", "T", "G", "0.500000"),
}

BUBBLES = [
    ("chr1", 1000, [
        ("A", "G", ["v1"]),
        ("A", "T", ["v1", "v2"]),   
        ("A", "C", ["v3"]),
    ]),
    ("chr1", 2000, [
        ("G", "A", ["v4"]),
        ("G", "C", ["v5"]),
        ("G", "T", ["v4", "v5"]),   
    ]),
    ("chr1", 3000, [
        ("T", "C", ["v6"]),
    ]),
    ("chr1", 4000, [
        ("C", "A", ["v7"]),
    ]),
    ("chr1", 5000, [
        ("A", "G", ["v8"]),
        ("A", "C", ["v9"]),
        ("A", "T", ["v8", "v9"]),
        ("A", "TT", ["v8"]),
    ]),
]

GT_CHOICES = ["1|0", "0|1", "1|1", "0|0"]

def sample_gp(rng):
    a = rng.randint(20, 960)
    b = rng.randint(20, 1000 - a - 20)
    c = 1000 - a - b
    return ("%.3f" % (a / 1000.0), "%.3f" % (b / 1000.0), "%.3f" % (c / 1000.0))

def build_group(bubble, rng):
    chrom, pos, alts = bubble
    group = []
    
    for (ref, alt, aids) in alts:
        samples = []
        for _ in range(NUM_SAMPLES):
            if pos == 3000:
                gt = "0|0"
                gp = ("1.000", "0.000", "0.000")
            elif pos == 4000:
                gt = "0|0"
                gp = ("0.9998", "0.0002", "0.0000")
            else:
                gt = rng.choice(GT_CHOICES)
                gp = sample_gp(rng)
                
            samples.append((gt, (float(gp[0]), float(gp[1]), float(gp[2]))))
            
        group.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "atomic_ids": aids, "samples": samples,
        })
    return group

MARGIN = 0.03   

def generate(seed=20240717):
    rng = random.Random(seed)
    final_groups = []
    # FIX: Isolate rejection sampling to individual bubbles
    for b in BUBBLES:
        while True:
            g = build_group(b, rng)
            ok = True
            for max_alleles in MAX_ALLELES_CASES:
                _, slack = process_group(g, ID_BUFFER, max_alleles)
                if slack < MARGIN:
                    ok = False
                    break
            if ok:
                final_groups.append(g)
                break
    return final_groups

CONTIG = "##contig=<ID=chr1,length=100000>"
SAMPLE_NAMES = ["S%d" % (i + 1) for i in range(NUM_SAMPLES)]

MAIN_HEADER = [
    "##fileformat=VCFv4.2",
    CONTIG,
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">',
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="genotype probs">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLE_NAMES),
]

EXPECTED_ADDED_HEADERS = [
    '##INFO=<ID=ID,Number=1,Type=String,Description="Atomic variant ID">',
    '##INFO=<ID=RAF,Number=A,Type=Float,Description="ALT allele frequency in the reference panel">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency in the target samples (mean dosage / 2)">',
    '##INFO=<ID=INFO,Number=A,Type=Float,Description="Atomic IMPUTE INFO recalculated from bubble-path GLIMPSE2 GP, defined as 1 when AF is 0 or 1">',
    '##INFO=<ID=N_PATHS,Number=1,Type=Integer,Description="Number of bubble paths that contain this variant">',
    '##INFO=<ID=N_PATHS_TOTAL,Number=1,Type=Integer,Description="Number of total non-reference paths in the bubble that contains this variant">'
]

def write_main(path, groups):
    with open(path, "w") as f:
        for h in MAIN_HEADER:
            f.write(h + "\n")
        for g in groups:
            for rec in g:
                cells = []
                for (gt, gp) in rec["samples"]:
                    cells.append("%s:%s,%s,%s" % (gt, "%.3f" % gp[0], "%.3f" % gp[1], "%.3f" % gp[2]))
                f.write("\t".join([rec["chrom"], str(rec["pos"]), ".", rec["ref"], rec["alt"], ".", ".", ".", "GT:GP"] + cells) + "\n")

def write_sites(path, groups):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n" + CONTIG + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for g in groups:
            for rec in g:
                ids = rec["atomic_ids"]
                idfield = "ID=" + ("," if len(ids) <= 1 else ":").join(ids)
                f.write("\t".join([rec["chrom"], str(rec["pos"]), ".", rec["ref"], rec["alt"], ".", ".", idfield]) + "\n")

def write_ids(path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n" + CONTIG + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key, (pos, orig, ref, alt, raf) in sorted(ID_BUFFER.items(), key=lambda kv: (kv[1][0], kv[0])):
            f.write("\t".join(["chr1", str(pos), orig, ref, alt, ".", ".", f"ID={key};AF={raf}"]) + "\n")

def expected_body(groups, max_alleles):
    out = []
    for g in groups:
        lines, _ = process_group(g, ID_BUFFER, max_alleles)
        out.extend(lines)
    return out

def write_expected(path, groups, max_alleles):
    with open(path, "w") as f:
        for h in MAIN_HEADER:
            if h.startswith("#CHROM"):
                for ah in EXPECTED_ADDED_HEADERS:
                    f.write(ah + "\n")
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

    for m in MAX_ALLELES_CASES:
        write_expected(os.path.join(exp, f"max{m}.txt"), groups, m)

if __name__ == "__main__":
    main()
