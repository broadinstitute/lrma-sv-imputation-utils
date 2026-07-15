#!/usr/bin/env python3
"""Generate deterministic fixtures and expected outputs for extract_bubble_PLs.

This is a *faithful re-implementation* of the tool's matching / PL logic (which
is entirely integer + string, so it is bit-exact against the Rust binary). It
writes:

    panel.vcf                 bi-allelic panel bubbles (INFO/ID carries the id)
    input.vcf                 gVCF-style input (multi-allelic PL, LPL, ref block)
    samples_S2.txt            sample-subset list
    expected/gvcf.txt         bcftools-query body, mode=gvcf, all samples
    expected/joint.txt        mode=joint, all samples
    expected/samples_S2.txt   mode=gvcf, --samples S2
    expected/region.txt       mode=gvcf, --region chr1:1400-2500

The generator hard-codes a hand-derived expectation for the gvcf case and
asserts the oracle reproduces it, so a mistake in either the fixtures or the
oracle is caught here rather than surfacing as a confusing CI diff.

Expected files are formatted to match:
    bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT:%PL]\\n'
"""

import math
import os
import sys

CAP_PL = 30
SCALE_PL = 5.0

# --------------------------------------------------------------------------- #
# Faithful oracle of the tool's arithmetic
# --------------------------------------------------------------------------- #
def process_pl(v: int) -> int:
    # ((v as f64) / scale) as i32  == truncation toward zero, then cap.
    scaled = math.trunc(v / SCALE_PL)
    return CAP_PL if scaled > CAP_PL else scaled


def minrep(pos, ref, alt):
    """Right-trim then left-trim shared bases; symbolic ALTs pass through."""
    if alt.startswith("<"):
        return (pos, ref, alt)
    r, a = list(ref), list(alt)
    while r and a and r[-1] == a[-1]:
        r.pop(); a.pop()
    while r and a and r[0] == a[0]:
        r.pop(0); a.pop(0); pos += 1
    return (pos, "".join(r), "".join(a))


def render_gt(a0, a1):
    def m(x):
        return "." if x is None else str(x)
    return f"{m(a0)}/{m(a1)}"   # always unphased in output


def emit_site(bubble, records, mode, samples):
    """Return list[(gt, pl)] per sample, or 'DROP', or None (skipped)."""
    if bubble["alt"] == ".":
        return None  # <2 alleles -> skipped entirely
    minp = minrep(bubble["pos"], bubble["ref"], bubble["alt"])

    # 1) variant match
    for rec in records:
        if rec["chrom"] != bubble["chrom"] or rec["ref_block"]:
            continue
        n = len(rec["alleles"])
        stride = n * (n + 1) // 2
        for k in range(1, n):
            alt = rec["alleles"][k]
            if alt == "<NON_REF>":
                continue
            mg = minrep(rec["pos"], rec["alleles"][0], alt)
            if mg[1] == "" and mg[2] == "":
                continue
            if mg == minp:
                idx0, idx1, idx2 = 0, k * (k + 1) // 2, k * (k + 1) // 2 + k
                out = []
                for s in samples:
                    a0, a1 = rec["gt"][s]
                    def mapa(x):
                        if x == 0: return 0
                        if x == k: return 1
                        return None
                    gt = render_gt(mapa(a0), mapa(a1))
                    src = rec["lpl"][s] if s in rec.get("lpl", {}) else rec["pl"][s]
                    if stride > idx2:
                        pl = [process_pl(src[idx0]), process_pl(src[idx1]),
                              process_pl(src[idx2])]
                    else:
                        pl = [0, 0, 0]
                    out.append((gt, ",".join(map(str, pl))))
                return out

    # 2) no match
    if mode == "joint":
        return "DROP"

    # gvcf: covering reference block -> hom-ref from GQ
    for rec in records:
        if rec["chrom"] != bubble["chrom"] or not rec["ref_block"]:
            continue
        site_end = minp[0] + len(minp[1]) - 1
        if rec["pos"] <= minp[0] and rec["end"] >= site_end:
            out = []
            for s in samples:
                gq = rec["gq"][s]
                out.append(("0/0", f"0,{process_pl(gq)},{process_pl(gq * 2)}"))
            return out

    # gvcf: uncovered -> missing
    return [("./.", "0,0,0") for _ in samples]


def build_expected(panel, records, mode, samples, region=None):
    lines = []
    for b in panel:
        if region is not None:
            chrom, lo, hi = region
            if b["chrom"] != chrom or b["pos"] < lo or b["pos"] > hi:
                continue
        res = emit_site(b, records, mode, samples)
        if res is None or res == "DROP":
            continue
        cols = [b["chrom"], str(b["pos"]), b["id"], b["ref"], b["alt"]]
        for gt, pl in res:
            cols.append(f"{gt}:{pl}")
        lines.append("\t".join(cols))
    return "\n".join(lines) + ("\n" if lines else "")


# --------------------------------------------------------------------------- #
# Fixtures (hand-designed to exercise the branches)
# --------------------------------------------------------------------------- #
SAMPLES = ["S1", "S2"]

PANEL = [
    {"chrom": "chr1", "pos": 1000, "ref": "A", "alt": "G", "id": "bub1"},  # multiallelic PL match
    {"chrom": "chr1", "pos": 1500, "ref": "C", "alt": "T", "id": "bub2"},  # LPL-preferred match
    {"chrom": "chr1", "pos": 2000, "ref": "G", "alt": "A", "id": "bub3"},  # ref-block hom-ref / joint-drop
    {"chrom": "chr1", "pos": 3000, "ref": "T", "alt": "C", "id": "bub4"},  # uncovered missing / joint-drop
    {"chrom": "chr1", "pos": 4000, "ref": "A", "alt": ".", "id": "bub5"},  # no ALT -> skipped
    {"chrom": "chr2", "pos": 1000, "ref": "A", "alt": "T", "id": "bub6"},  # different contig match
]

# Input records. gt[sample] = (allele0, allele1) as integer indices.
RECORDS = [
    {   # multi-allelic variant, standard PL (stride 6)
        "chrom": "chr1", "pos": 1000, "ref_block": False,
        "alleles": ["A", "G", "C"],
        "gt": {"S1": (0, 1), "S2": (1, 2)},
        "pl": {"S1": [0, 10, 50, 20, 60, 70], "S2": [30, 20, 0, 45, 5, 60]},
    },
    {   # bi-allelic variant carrying BOTH LPL and PL; LPL must win
        "chrom": "chr1", "pos": 1500, "ref_block": False,
        "alleles": ["C", "T"],
        "gt": {"S1": (0, 0), "S2": (0, 1)},
        "lpl": {"S1": [0, 15, 99], "S2": [20, 0, 40]},
        "pl":  {"S1": [5, 5, 5],   "S2": [5, 5, 5]},   # decoy: would give 1,1,1
    },
    {   # reference block covering chr1:2000 only
        "chrom": "chr1", "pos": 1900, "ref_block": True, "end": 2100,
        "alleles": ["A", "<NON_REF>"],
        "gq": {"S1": 35, "S2": 12},
    },
    {   # chr2 variant
        "chrom": "chr2", "pos": 1000, "ref_block": False,
        "alleles": ["A", "T"],
        "gt": {"S1": (0, 1), "S2": (0, 0)},
        "pl": {"S1": [0, 8, 88], "S2": [0, 50, 120]},
    },
]

# Hand-derived expectation for the gvcf/all-samples case (double-entry check).
HAND_GVCF = (
    "chr1\t1000\tbub1\tA\tG\t0/1:0,2,10\t1/.:6,4,0\n"
    "chr1\t1500\tbub2\tC\tT\t0/0:0,3,19\t0/1:4,0,8\n"
    "chr1\t2000\tbub3\tG\tA\t0/0:0,7,14\t0/0:0,2,4\n"
    "chr1\t3000\tbub4\tT\tC\t./.:0,0,0\t./.:0,0,0\n"
    "chr2\t1000\tbub6\tA\tT\t0/1:0,1,17\t0/0:0,10,24\n"
)


# --------------------------------------------------------------------------- #
# VCF writers
# --------------------------------------------------------------------------- #
CONTIGS = [("chr1", 100000), ("chr2", 100000)]


def write_panel(path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c, ln in CONTIGS:
            f.write(f"##contig=<ID={c},length={ln}>\n")
        f.write('##INFO=<ID=ID,Number=1,Type=String,Description="Panel bubble ID">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for b in PANEL:
            f.write("\t".join([b["chrom"], str(b["pos"]), ".", b["ref"], b["alt"],
                               ".", ".", f"ID={b['id']}"]) + "\n")


def write_input(path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c, ln in CONTIGS:
            f.write(f"##contig=<ID={c},length={ln}>\n")
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        f.write('##ALT=<ID=NON_REF,Description="Non-ref block">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">\n')
        f.write('##FORMAT=<ID=LPL,Number=.,Type=Integer,Description="Local PL">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(SAMPLES) + "\n")
        # Emit records sorted by (chrom, pos) as a valid VCF requires.
        rows = []
        for r in RECORDS:
            if r["ref_block"]:
                info = f"END={r['end']}"
                fmt = "GT:GQ"
                cells = []
                for s in SAMPLES:
                    cells.append(f"0/0:{r['gq'][s]}")
                alt = r["alleles"][1]
            else:
                info = "."
                has_lpl = "lpl" in r
                fmt = "GT:PL" + (":LPL" if has_lpl else "")
                cells = []
                for s in SAMPLES:
                    gt = "/".join(str(x) for x in r["gt"][s])
                    pl = ",".join(map(str, r["pl"][s]))
                    cell = f"{gt}:{pl}"
                    if has_lpl:
                        cell += ":" + ",".join(map(str, r["lpl"][s]))
                    cells.append(cell)
                alt = ",".join(r["alleles"][1:])
            rows.append((r["chrom"], r["pos"],
                         "\t".join([r["chrom"], str(r["pos"]), ".", r["alleles"][0],
                                    alt, ".", ".", info, fmt] + cells)))
        for _, _, line in sorted(rows, key=lambda t: (t[0], t[1])):
            f.write(line + "\n")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    exp = os.path.join(here, "expected")
    os.makedirs(exp, exist_ok=True)

    write_panel(os.path.join(here, "panel.vcf"))
    write_input(os.path.join(here, "input.vcf"))
    with open(os.path.join(here, "samples_S2.txt"), "w") as f:
        f.write("S2\n")

    gvcf = build_expected(PANEL, RECORDS, "gvcf", SAMPLES)
    if gvcf != HAND_GVCF:
        sys.stderr.write("ORACLE SELF-CHECK FAILED (gvcf)\n--- oracle ---\n"
                         + gvcf + "\n--- hand ---\n" + HAND_GVCF + "\n")
        sys.exit(2)

    with open(os.path.join(exp, "gvcf.txt"), "w") as f:
        f.write(gvcf)
    with open(os.path.join(exp, "joint.txt"), "w") as f:
        f.write(build_expected(PANEL, RECORDS, "joint", SAMPLES))
    with open(os.path.join(exp, "samples_S2.txt"), "w") as f:
        f.write(build_expected(PANEL, RECORDS, "gvcf", ["S2"]))
    with open(os.path.join(exp, "region.txt"), "w") as f:
        f.write(build_expected(PANEL, RECORDS, "gvcf", SAMPLES,
                               region=("chr1", 1400, 2500)))

    print("extract-bubble-PLs fixtures + expected written; oracle self-check OK")


if __name__ == "__main__":
    main()
