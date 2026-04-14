#!/usr/bin/env bash

# Simple helper script to generate visualization-ready results
# for the P12-noXY sample and hook them into the Witness UI.
#
# Usage:
#   1) Make it executable:
#        chmod +x PhyloWGS_runtime_analysis/visualize_P12-noXY.sh
#   2) Run:
#        ./PhyloWGS_runtime_analysis/visualize_P12-noXY.sh
#   3) Then start a web server from the project root, e.g.:
#        cd /path/to/PhyloWGS
#        python3 -m http.server 8000
#      and open:
#        http://localhost:8000/witness/index.html

set -euo pipefail

# Resolve project root (this script lives in PhyloWGS_runtime_analysis/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RESULT_DIR="${PROJECT_ROOT}/PhyloWGS_runtime_analysis/results_phylowgs/P12-noXY"
WITNESS_DATA_RUN="dream"
WITNESS_DATA_DIR="${PROJECT_ROOT}/witness/data/${WITNESS_DATA_RUN}"

PYTHON2_BIN="${PYTHON2_BIN:-python2}"

echo "Project root     : ${PROJECT_ROOT}"
echo "Result directory : ${RESULT_DIR}"
echo "Witness data dir : ${WITNESS_DATA_DIR}"
echo

if ! command -v "${PYTHON2_BIN}" >/dev/null 2>&1; then
  echo "ERROR: python2 not found in PATH (or via \$PYTHON2_BIN)." >&2
  echo "       Please install Python 2 or point PYTHON2_BIN to a Python 2 executable." >&2
  exit 1
fi

TREES_ZIP="${RESULT_DIR}/trees.zip"
if [[ ! -f "${TREES_ZIP}" ]]; then
  echo "ERROR: Cannot find trees.zip at ${TREES_ZIP}." >&2
  echo "       Make sure PhyloWGS has completed for P12-noXY." >&2
  exit 1
fi

mkdir -p "${WITNESS_DATA_DIR}"

TREES_SUMM_OUT="${WITNESS_DATA_DIR}/P12-noXY.summ.json.gz"
MUTLIST_OUT="${WITNESS_DATA_DIR}/P12-noXY.muts.json.gz"
MUTASS_OUT="${WITNESS_DATA_DIR}/P12-noXY.mutass.zip"

echo "Generating JSON / mutass results for P12-noXY..."
echo "  Summary : ${TREES_SUMM_OUT}"
echo "  Mutlist : ${MUTLIST_OUT}"
echo "  Mutass  : ${MUTASS_OUT}"
echo

(
  cd "${PROJECT_ROOT}"
  "${PYTHON2_BIN}" write_results.py \
    P12-noXY \
    "${TREES_ZIP}" \
    "${TREES_SUMM_OUT}" \
    "${MUTLIST_OUT}" \
    "${MUTASS_OUT}"
)

echo
echo "Decompressing for Witness UI (browser needs plain JSON)..."
SUMM_JSON="${WITNESS_DATA_DIR}/P12-noXY.summ.json"
MUTS_JSON="${WITNESS_DATA_DIR}/P12-noXY.muts.json"
MUTASS_DIR="${WITNESS_DATA_DIR}/P12-noXY.mutass"

gunzip -c "${TREES_SUMM_OUT}" > "${SUMM_JSON}"
gunzip -c "${MUTLIST_OUT}" > "${MUTS_JSON}"

rm -rf "${MUTASS_DIR}"
mkdir -p "${MUTASS_DIR}"
unzip -o -j "${MUTASS_OUT}" -d "${MUTASS_DIR}" >/dev/null 2>&1

echo "  Created: ${SUMM_JSON}"
echo "  Created: ${MUTS_JSON}"
echo "  Created: ${MUTASS_DIR}/ (from mutass.zip)"
echo
echo "Updating Witness index (witness/data/index.json)..."
(
  cd "${PROJECT_ROOT}/witness"
  "${PYTHON2_BIN}" index_data.py
)

echo
echo "Done."
echo
echo "Next steps to view the results:"
echo "  1) From the project root (${PROJECT_ROOT}), start a simple web server, e.g.:"
echo "       python3 -m http.server 8000"
echo "  2) Open this URL in your browser:"
echo "       http://localhost:8000/witness/index.html"
echo "  3) In the Witness UI, look for run \"${WITNESS_DATA_RUN}\" and dataset \"P12-noXY\""
echo "     to inspect trees, clusters and VAF plots."

