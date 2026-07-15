#!/usr/bin/env bash
# Run every per-tool integration suite and aggregate the result.
#
# Usage:
#   resources/tests/run_all.sh [--bin-dir DIR]
#
# Binary discovery order (per tool, see lib/harness.sh):
#   1. tool-specific env override (EXTRACT_BUBBLE_PLS_BIN / POP_GLIMPSE2_BIN /
#      PASTE_VCFS_BIN),
#   2. $BIN_DIR (this flag, or the env var; default target/release),
#   3. $PATH.
#
# Inside the final Docker image the binaries live in /usr/local/bin, so CI runs
#   resources/tests/run_all.sh --bin-dir /usr/local/bin
# with the tests bind-mounted (never baked into the image).

set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin-dir) export BIN_DIR="$2"; shift 2 ;;
    --bin-dir=*) export BIN_DIR="${1#*=}"; shift ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

SUITES=(
  "extract-bubble-PLs/run.sh"
  "pop-glimpse2/run.sh"
  "paste-vcfs/run.sh"
)

overall=0
for s in "${SUITES[@]}"; do
  echo
  echo "======================================================================"
  echo "  suite: $s"
  echo "======================================================================"
  if bash "${HERE}/${s}"; then
    echo ">> ${s}: OK"
  else
    echo ">> ${s}: FAILED"
    overall=1
  fi
done

echo
if [[ $overall -eq 0 ]]; then
  echo "ALL SUITES PASSED"
else
  echo "ONE OR MORE SUITES FAILED"
fi
exit $overall
