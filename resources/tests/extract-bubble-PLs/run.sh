#!/usr/bin/env bash
# Integration tests for extract_bubble_PLs.
#
# Strategy: the arithmetic of this tool is entirely integer + string, so a
# faithful Python re-implementation (generate.py) is a bit-exact oracle. We
# generate fixtures + expected once, build the indexed BCF/VCF.GZ inputs the
# binary needs, run the binary, and compare `bcftools query` bodies exactly.

set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/harness.sh
source "${HERE}/../lib/harness.sh"

require_tools
EXTRACT="$(bin_extract)"
info "using binary: $EXTRACT"

# 1) (Re)generate deterministic fixtures + expected outputs from the oracle.
python3 "${HERE}/generate.py"

SCRATCH="$(new_scratch)"

# 2) Build the indexed inputs the tool opens via IndexedReader::from_path.
PANEL_BCF="${SCRATCH}/panel.bcf"
INPUT_GZ="${SCRATCH}/input.vcf.gz"
make_bcf   "${HERE}/panel.vcf" "$PANEL_BCF"
make_vcfgz "${HERE}/input.vcf" "$INPUT_GZ"

# Query format must mirror generate.py's expected files exactly.
QFMT='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%PL]\n'

run_case() {
  # run_case <name> <expected-file> <mode> [extra args...]
  local name="$1" expected="$2" mode="$3"; shift 3
  local out="${SCRATCH}/${name}.bcf"
  "$EXTRACT" "$mode" "$PANEL_BCF" "$INPUT_GZ" "$out" "$@" \
    >"${SCRATCH}/${name}.log" 2>&1 || {
      log "${C_RED}FAIL${C_RST} ${name}: binary exited non-zero"
      sed 's/^/      /' "${SCRATCH}/${name}.log" >&2
      TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1)); return 0
    }
  query_body "$out" "$QFMT" > "${SCRATCH}/${name}.actual"
  assert_files_equal "$name" "${SCRATCH}/${name}.actual" "$expected" || true
}

# 3) Cases — each exercises a distinct branch of the matcher.
#    gvcf, all samples: multiallelic alt-k pick, LPL-over-PL, ref-block hom-ref,
#                       uncovered ./., <2-allele skip, second contig.
run_case gvcf        "${HERE}/expected/gvcf.txt"        gvcf
#    joint, all samples: matched rows emitted, unmatched rows DROPPED.
run_case joint       "${HERE}/expected/joint.txt"       joint
#    --samples subset restricts the sample columns.
run_case samples_S2  "${HERE}/expected/samples_S2.txt"  gvcf --samples "${HERE}/samples_S2.txt"
#    --region clips the panel iteration to a 0-based-guarded window.
run_case region      "${HERE}/expected/region.txt"      gvcf --region chr1:1400-2500

finish
