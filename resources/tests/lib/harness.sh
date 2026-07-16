#!/usr/bin/env bash
# Shared helpers for the Rust-tool integration tests.
#
# Responsibilities:
#   * locate the three compiled binaries (env override -> BIN_DIR -> PATH),
#   * turn VCF text fixtures into indexed BCF / VCF.GZ at run time,
#   * run bcftools and normalise output for exact-match comparison,
#   * provide pass/fail assertions and a summary that sets the exit code.
#
# The tests never bake binaries or data into anything; fixtures are generated
# into a scratch dir and removed on exit.

set -Eeuo pipefail

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
# REPO_ROOT is two levels up from resources/tests/lib.
HARNESS_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HARNESS_LIB_DIR}/../../.." && pwd)"
# Where cargo puts release binaries for the workspace.
BIN_DIR="${BIN_DIR:-${REPO_ROOT}/target/release}"

# --------------------------------------------------------------------------
# Counters
# --------------------------------------------------------------------------
TESTS_RUN=0
TESTS_FAILED=0
TESTS_SKIPPED=0

C_RED=$'\033[31m'; C_GRN=$'\033[32m'; C_YEL=$'\033[33m'; C_DIM=$'\033[2m'; C_RST=$'\033[0m'
if [[ ! -t 1 ]]; then C_RED=; C_GRN=; C_YEL=; C_DIM=; C_RST=; fi

log()  { printf '%s\n' "$*" >&2; }
info() { printf '%s[info]%s %s\n' "$C_DIM" "$C_RST" "$*" >&2; }

# --------------------------------------------------------------------------
# Tool dependency checks
# --------------------------------------------------------------------------
require_tools() {
  local missing=0 t
  # Deliberately bcftools-only: the final runtime image ships bcftools + python3
  # and nothing else, so every bgzip/index operation goes through bcftools
  # (view -Oz / index / index -t) rather than standalone bgzip / tabix.
  for t in bcftools python3 diff; do
    if ! command -v "$t" >/dev/null 2>&1; then
      log "${C_RED}[fatal]${C_RST} required tool not found on PATH: $t"
      missing=1
    fi
  done
  [[ $missing -eq 0 ]] || exit 3
}

# Resolve a binary from a list of candidate names, honouring an override var.
# usage: resolve_bin OVERRIDE_VAR "cand1 cand2 ..."
resolve_bin() {
  local override_var="$1" candidates="$2" c p
  local override="${!override_var:-}"
  if [[ -n "$override" ]]; then
    if [[ -x "$override" ]]; then printf '%s\n' "$override"; return 0; fi
    log "${C_RED}[fatal]${C_RST} \$$override_var=$override is not executable"; exit 3
  fi
  for c in $candidates; do
    if [[ -x "$BIN_DIR/$c" ]]; then printf '%s\n' "$BIN_DIR/$c"; return 0; fi
  done
  for c in $candidates; do
    if p="$(command -v "$c" 2>/dev/null)"; then printf '%s\n' "$p"; return 0; fi
  done
  log "${C_RED}[fatal]${C_RST} could not find any of: $candidates"
  log "        looked in BIN_DIR=$BIN_DIR and on PATH."
  log "        set \$$override_var to point at the binary."
  exit 3
}

bin_extract()     { resolve_bin EXTRACT_BUBBLE_PLS_BIN "extract_bubble_PLs extract-bubble-PLs"; }
bin_popglimpse()  { resolve_bin POP_GLIMPSE2_BIN       "pop-glimpse2"; }
bin_paste()       { resolve_bin PASTE_VCFS_BIN         "paste-vcfs paste_vcfs"; }

# --------------------------------------------------------------------------
# Fixture builders
# --------------------------------------------------------------------------
# VCF text -> compressed+CSI-indexed BCF.
make_bcf() {
  local in_vcf="$1" out_bcf="$2"
  bcftools view -Ob -o "$out_bcf" "$in_vcf" >/dev/null 2>&1
  bcftools index -f "$out_bcf"            # writes .csi
}

# VCF text -> bgzipped VCF with both TBI and CSI (IndexedReader auto-detects).
# Uses bcftools for the compression and indexing so no standalone bgzip/tabix
# binary is required (they are absent from the runtime image).
make_vcfgz() {
  local in_vcf="$1" out_gz="$2"
  bcftools view -Oz -o "$out_gz" "$in_vcf" >/dev/null 2>&1
  bcftools index -f -t "$out_gz"          # writes .tbi
  bcftools index -f -c "$out_gz"          # writes .csi (some readers prefer it)
}

# --------------------------------------------------------------------------
# Normalised readouts
# --------------------------------------------------------------------------
# Deterministic body of a BCF/VCF for exact-match comparison.
query_body() {
  local file="$1" fmt="$2"
  bcftools query -f "$fmt" "$file"
}

# Header lines that survive across bcftools versions (drops volatile provenance).
norm_header() {
  local file="$1"
  bcftools view -h "$file" \
    | grep -Ev '^##(bcftools|fileDate|source|GATKCommandLine|reference=)' || true
}

# --------------------------------------------------------------------------
# Assertions
# --------------------------------------------------------------------------
# Compare "actual" text (arg 2, a file) to an expected file (arg 3).
assert_files_equal() {
  local name="$1" actual="$2" expected="$3"
  TESTS_RUN=$((TESTS_RUN + 1))
  if [[ ! -f "$expected" ]]; then
    log "${C_RED}FAIL${C_RST} ${name}: expected file missing: $expected"
    TESTS_FAILED=$((TESTS_FAILED + 1)); return 1
  fi
  if diff -u "$expected" "$actual" > "${actual}.diff" 2>&1; then
    log "${C_GRN}PASS${C_RST} ${name}"
    rm -f "${actual}.diff"
    return 0
  else
    log "${C_RED}FAIL${C_RST} ${name}: output differs from expected"
    sed 's/^/      /' "${actual}.diff" >&2
    TESTS_FAILED=$((TESTS_FAILED + 1)); return 1
  fi
}

# Assert a command fails (non-zero) and its stderr matches a regex.
assert_fails_matching() {
  local name="$1" pattern="$2"; shift 2
  TESTS_RUN=$((TESTS_RUN + 1))
  local err rc=0
  err="$("$@" 2>&1 1>/dev/null)" || rc=$?
  if [[ $rc -eq 0 ]]; then
    log "${C_RED}FAIL${C_RST} ${name}: command unexpectedly succeeded"
    TESTS_FAILED=$((TESTS_FAILED + 1)); return 1
  fi
  if grep -Eq "$pattern" <<<"$err"; then
    log "${C_GRN}PASS${C_RST} ${name} (rejected as expected)"
    return 0
  fi
  log "${C_RED}FAIL${C_RST} ${name}: exited $rc but stderr did not match /$pattern/"
  sed 's/^/      /' <<<"$err" >&2
  TESTS_FAILED=$((TESTS_FAILED + 1)); return 1
}

assert_contains() {
  local name="$1" file="$2" pattern="$3"
  TESTS_RUN=$((TESTS_RUN + 1))
  if grep -Eq "$pattern" "$file"; then
    log "${C_GRN}PASS${C_RST} ${name}"
    return 0
  fi
  log "${C_RED}FAIL${C_RST} ${name}: /$pattern/ not found in $file"
  TESTS_FAILED=$((TESTS_FAILED + 1)); return 1
}

skip() {
  local name="$1" why="$2"
  TESTS_SKIPPED=$((TESTS_SKIPPED + 1))
  log "${C_YEL}SKIP${C_RST} ${name}: ${why}"
}

finish() {
  log ""
  log "  ran ${TESTS_RUN}  |  failed ${TESTS_FAILED}  |  skipped ${TESTS_SKIPPED}"
  [[ $TESTS_FAILED -eq 0 ]] || exit 1
  exit 0
}

# Scratch space. IMPORTANT: the root is created here, at source time, in the
# suite's main shell -- NOT inside new_scratch(). new_scratch() is called as
# `SCRATCH="$(new_scratch)"`, i.e. in a command-substitution subshell; a
# `trap ... EXIT` registered inside that function would fire the moment the
# subshell exits and delete the directory before the caller could use it. By
# owning the root (and its single EXIT trap) in the main shell, subdirectories
# handed out by new_scratch() survive until the suite process itself exits.
_SCRATCH_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/rusttest.XXXXXX")"
# shellcheck disable=SC2064
trap 'rm -rf "${_SCRATCH_ROOT:-}"' EXIT

# Hand out a fresh subdirectory under the persistent root.
new_scratch() {
  local d; d="$(mktemp -d "${_SCRATCH_ROOT}/s.XXXXXX")"
  printf '%s\n' "$d"
}

# Choose a WRITABLE directory for a suite's generated fixtures/expected. Returns
# the given source directory if it is writable (native `make test` against a
# checkout), otherwise a fresh scratch subdirectory (e.g. `make docker-test`,
# where resources/tests is bind-mounted read-only). The generators honour this
# via the GEN_OUT_DIR environment variable.
writable_datadir() {
  local src="$1" probe
  probe="${src}/.wtest.$$"
  if ( : > "$probe" ) 2>/dev/null; then
    rm -f "$probe"
    printf '%s\n' "$src"
  else
    mktemp -d "${_SCRATCH_ROOT}/data.XXXXXX"
  fi
}
