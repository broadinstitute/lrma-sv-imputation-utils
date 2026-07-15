#!/usr/bin/env bash
# Integration tests for pop-glimpse2-joint-opt.
#
# The tool reads a multiallelic VCF on STDIN plus two file arguments
# (biallelic-ID VCF, sites VCF) and writes the projected VCF to STDOUT. No
# bcftools / indexing is involved -- it is plain text (or .gz) throughout.
#
# generate.py is a faithful f32 re-implementation whose fixtures are chosen by
# rejection sampling to sit far from every rounding / threshold boundary, so
# the Rust f32 output is byte-identical to the expected files.

set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/harness.sh
source "${HERE}/../lib/harness.sh"

require_tools
POPG="$(bin_popglimpse)"
info "using binary: $POPG"

python3 "${HERE}/generate.py"

SCRATCH="$(new_scratch)"

run_case() {
  # run_case <name> <expected-file> [max_alleles arg...]
  local name="$1" expected="$2"; shift 2
  local actual="${SCRATCH}/${name}.actual"
  if "$POPG" "${HERE}/ids.vcf" "${HERE}/sites.vcf" "$@" \
        < "${HERE}/main.vcf" > "$actual" 2>"${SCRATCH}/${name}.log"; then
    assert_files_equal "$name" "$actual" "$expected" || true
  else
    log "${C_RED}FAIL${C_RST} ${name}: binary exited non-zero"
    sed 's/^/      /' "${SCRATCH}/${name}.log" >&2
    TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1))
  fi
}

# Default max_alleles (10) >= alleles-per-bubble: full projection, every allele
# contributes to the odds-normalisation.
run_case max10 "${HERE}/expected/max10.txt"
# max_alleles = 2 < 3 alleles: exercises the stable top-k score cut.
run_case max2  "${HERE}/expected/max2.txt" 2

# Lockstep guard: a sites file whose ALT disagrees with the main stream must
# abort with the synchronisation error rather than emit anything.
BAD_SITES="${SCRATCH}/bad_sites.vcf"
sed '0,/A\tG/s//A\tG_BAD/' "${HERE}/sites.vcf" > "$BAD_SITES" 2>/dev/null || cp "${HERE}/sites.vcf" "$BAD_SITES"
# Force a guaranteed mismatch on the first data ALT.
awk 'BEGIN{OFS="\t"} /^#/{print;next} {if(!done){$5=$5"X";done=1} print}' \
    "${HERE}/sites.vcf" > "$BAD_SITES"
assert_fails_matching "lockstep_mismatch_aborts" "Lockstep|do not match|synchron" \
    bash -c '"$1" "$2" "$3" < "$4"' _ \
        "$POPG" "${HERE}/ids.vcf" "$BAD_SITES" "${HERE}/main.vcf"

finish
