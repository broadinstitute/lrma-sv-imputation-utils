#!/usr/bin/env bash
# Integration tests for paste-vcfs.
#
# There is no arithmetic to oracle, so the expected output is *derived from the
# inputs*: we re-query each input with the same bcftools format string the
# comparison uses and horizontally join the per-sample blocks. Because both
# sides of the diff go through identical bcftools formatting, any float
# rendering is cancelled out and the match is exact.

set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/harness.sh
source "${HERE}/../lib/harness.sh"

require_tools
PASTE="$(bin_paste)"
info "using binary: $PASTE"

python3 "${HERE}/generate.py"

SCRATCH="$(new_scratch)"

# Build indexed BCFs (make_bcf writes .csi, needed by region mode).
for n in in0 in1 in2; do
  make_bcf "${HERE}/${n}.vcf" "${SCRATCH}/${n}.bcf"
done

# Comparison readouts. Leading columns come from the BASE (in0); sample blocks
# are each input's per-sample GT:PL:DS, tab-prefixed, joined in file order.
LEAD_FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/AF\n'
SAMP_FMT='[\t%GT:%PL:%DS]\n'

# make_expected <out> <region-or-empty> <base.bcf> [side.bcf ...]
make_expected() {
  local out="$1" region="$2"; shift 2
  local base="$1"; shift
  local ropt=()
  [[ -n "$region" ]] && ropt=(-r "$region")
  local blocks=() f
  blocks+=(<(bcftools query "${ropt[@]}" -f "$LEAD_FMT" "$base"))
  blocks+=(<(bcftools query "${ropt[@]}" -f "$SAMP_FMT" "$base"))
  for f in "$@"; do
    blocks+=(<(bcftools query "${ropt[@]}" -f "$SAMP_FMT" "$f"))
  done
  paste -d '' "${blocks[@]}" > "$out"
}

# run_case <name> <tool-args...> -- <inputs...>
#   tool-args : everything passed to paste-vcfs before -o/inputs (may include -r)
#   inputs    : the input BCFs, base first
run_case() {
  local name="$1"; shift
  local args=()
  while [[ "$1" != "--" ]]; do args+=("$1"); shift; done
  shift                     # drop the "--"
  local inputs=("$@")
  local out="${SCRATCH}/${name}.bcf"

  if "$PASTE" "${args[@]}" -o "$out" "${inputs[@]}" >"${SCRATCH}/${name}.log" 2>&1; then
    # recover an optional region from the tool args (for the expected build)
    local region="" i
    for ((i = 0; i < ${#args[@]}; i++)); do
      if [[ "${args[$i]}" == "-r" || "${args[$i]}" == "--region" ]]; then
        region="${args[$((i + 1))]}"
      fi
    done
    make_expected "${SCRATCH}/${name}.exp" "$region" "${inputs[@]}"
    bcftools query ${region:+-r "$region"} \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/AF[\t%GT:%PL:%DS]\n' \
        "$out" > "${SCRATCH}/${name}.actual"
    assert_files_equal "$name" "${SCRATCH}/${name}.actual" "${SCRATCH}/${name}.exp" || true
  else
    log "${C_RED}FAIL${C_RST} ${name}: binary exited non-zero"
    sed 's/^/      /' "${SCRATCH}/${name}.log" >&2
    TESTS_RUN=$((TESTS_RUN + 1)); TESTS_FAILED=$((TESTS_FAILED + 1))
  fi
}

# 2-file whole paste: A1,A2 + B1, base INFO/AF carried.
run_case paste2 --info AF --format GT,PL,DS -- "${SCRATCH}/in0.bcf" "${SCRATCH}/in1.bcf"
# 3-file whole paste: A1,A2 + B1 + C1,C2.
run_case paste3 --info AF --format GT,PL,DS -- "${SCRATCH}/in0.bcf" "${SCRATCH}/in1.bcf" "${SCRATCH}/in2.bcf"
# Region seek (needs .csi): restrict to chr1, dropping the chr2 site.
run_case region_chr1 -r chr1 --info AF --format GT,PL,DS -- "${SCRATCH}/in0.bcf" "${SCRATCH}/in1.bcf" "${SCRATCH}/in2.bcf"
# Bounded region on chr1.
run_case region_bounded -r chr1:200-500 --info AF --format GT,PL,DS -- "${SCRATCH}/in0.bcf" "${SCRATCH}/in1.bcf"

# --- panic paths -------------------------------------------------------------
# Site mismatch: a side file whose ALT differs must abort.
awk 'BEGIN{OFS="\t"} /^#/{print;next} {if(!d){$5=$5"X";d=1} print}' \
    "${HERE}/in1.vcf" > "${SCRATCH}/in1_bad.vcf"
make_bcf "${SCRATCH}/in1_bad.vcf" "${SCRATCH}/in1_bad.bcf"
assert_fails_matching "site_mismatch_aborts" "REF/ALT mismatch|Position mismatch|ID mismatch" \
    "$PASTE" --format GT,PL,DS -o "${SCRATCH}/none1.bcf" "${SCRATCH}/in0.bcf" "${SCRATCH}/in1_bad.bcf"

# FORMAT tag present in base but absent in a side file must abort.
awk 'BEGIN{OFS="\t"} /^#/{
        if($0 ~ /ID=DS/) next; print; next
     }
     { n=split($0,a,"\t"); line="";
       for(i=1;i<=n;i++){ c=a[i];
         if(i==9){ c="GT:PL" }
         else if(i>9){ split(c,p,":"); c=p[1]":"p[2] }
         line=line (i>1?"\t":"") c }
       print line }' \
    "${HERE}/in1.vcf" > "${SCRATCH}/in1_nods.vcf"
make_bcf "${SCRATCH}/in1_nods.vcf" "${SCRATCH}/in1_nods.bcf"
assert_fails_matching "format_absent_aborts" "FORMAT tag|mismatch" \
    "$PASTE" --format GT,PL,DS -o "${SCRATCH}/none2.bcf" "${SCRATCH}/in0.bcf" "${SCRATCH}/in1_nods.bcf"

finish
