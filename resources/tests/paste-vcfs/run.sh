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
#
# Builds the expected pasted output by re-querying each input and horizontally
# joining the per-sample blocks. NOTE: we materialise each query into a temp
# file and then `paste` the files, rather than `paste`-ing bash process
# substitutions. Storing `<(...)` in an array does NOT work: bash sets up each
# /dev/fd entry but tears it down as soon as the producing command exits (which
# is immediate for these tiny inputs), so by the time `paste` runs the FDs are
# already gone ("paste: /dev/fd/63: No such file or directory").
make_expected() {
  local out="$1" region="$2"; shift 2
  local base="$1"; shift
  local ropt=()
  [[ -n "$region" ]] && ropt=(-r "$region")

  local td; td="$(mktemp -d "${SCRATCH}/exp.XXXXXX")"
  local files=() i=0 f

  # leading columns + INFO/AF come from the base
  bcftools query "${ropt[@]}" -f "$LEAD_FMT" "$base" > "${td}/00_lead"
  files+=("${td}/00_lead")
  # base's own sample block
  bcftools query "${ropt[@]}" -f "$SAMP_FMT" "$base" > "${td}/01_base"
  files+=("${td}/01_base")
  # each side file's sample block, in order
  for f in "$@"; do
    i=$((i + 1))
    local bf; bf="${td}/$(printf '%02d' $((i + 1)))_side"
    bcftools query "${ropt[@]}" -f "$SAMP_FMT" "$f" > "$bf"
    files+=("$bf")
  done

  paste -d '' "${files[@]}" > "$out"
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
    # The tool has already restricted its output to $region, so read the whole
    # output (no -r; it isn't indexed). Comparing that against the region-
    # filtered inputs also confirms the tool applied the restriction correctly.
    bcftools query \
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
