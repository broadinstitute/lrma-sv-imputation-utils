#!/usr/bin/env python3
"""Generate fixtures for paste-vcfs.

paste-vcfs concatenates sample columns across inputs that share *identical*
sites (same CHROM/POS/ID/REF/ALT, in the same order). There is no arithmetic to
oracle: the expected output is derived from the inputs themselves at run time
(re-query each input with the same bcftools format string and horizontally join
the per-sample blocks). This generator therefore only needs to lay down inputs
with shared sites, distinct sample names, and per-sample FORMAT values that are
easy to eyeball in a diff.

Inputs:
  in0.vcf  samples A1,A2   (BASE: its INFO/ID/QUAL/AF are the ones carried)
  in1.vcf  sample  B1
  in2.vcf  samples C1,C2

FORMAT tags: GT (Number=1), PL (Integer, Number=G), DS (Float, Number=1).
INFO tag:    AF (Float, Number=A) -- present in all, but only base's is kept.

Two contigs (chr1, chr2) so the region tests can restrict to chr1.
"""

import os

CONTIGS = [("chr1", 1000000), ("chr2", 1000000)]

# Shared sites: (chrom, pos, id, ref, alt, af)
SITES = [
    ("chr1", 100, "rs1", "A", "G", "0.10"),
    ("chr1", 250, "rs2", "C", "T", "0.25"),
    ("chr1", 400, "rs3", "G", "A", "0.40"),
    ("chr1", 900, "rs4", "T", "C", "0.05"),
    ("chr2", 150, "rs5", "A", "C", "0.33"),
]

HEADER_INFO_FORMAT = [
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">',
    '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">',
]


def cell(gt, pl, ds):
    return "%s:%s:%s" % (gt, ",".join(str(x) for x in pl), ds)


def sample_values(tag_seed, site_idx, sample_idx):
    """Deterministic, distinct-looking GT:PL:DS per (file,sample,site)."""
    gts = ["0/0", "0/1", "1/1", "0/1", "1/1"]
    gt = gts[(site_idx + sample_idx + tag_seed) % len(gts)]
    base = (tag_seed * 100) + (sample_idx * 10) + site_idx
    pl = [base % 90, (base + 13) % 90, (base + 27) % 90]
    ds = "%.2f" % (((base % 200) / 100.0))
    return cell(gt, pl, ds)


def write_vcf(path, samples, tag_seed, af_offset=0.0):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c, ln in CONTIGS:
            f.write("##contig=<ID=%s,length=%d>\n" % (c, ln))
        for h in HEADER_INFO_FORMAT:
            f.write(h + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n")
        for si, (chrom, pos, vid, ref, alt, af) in enumerate(SITES):
            af_val = "%.2f" % (float(af) + af_offset)
            cells = [sample_values(tag_seed, si, sj) for sj in range(len(samples))]
            f.write("\t".join([chrom, str(pos), vid, ref, alt, ".", ".",
                               "AF=%s" % af_val, "GT:PL:DS"] + cells) + "\n")


def main():
    here = os.environ.get("GEN_OUT_DIR") or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(here, exist_ok=True)
    write_vcf(os.path.join(here, "in0.vcf"), ["A1", "A2"], tag_seed=0, af_offset=0.0)
    write_vcf(os.path.join(here, "in1.vcf"), ["B1"],       tag_seed=1, af_offset=0.5)
    write_vcf(os.path.join(here, "in2.vcf"), ["C1", "C2"], tag_seed=2, af_offset=0.3)
    print("paste-vcfs fixtures written (in0/in1/in2; shared sites, distinct samples)")


if __name__ == "__main__":
    main()
