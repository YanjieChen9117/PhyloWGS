#!/usr/bin/env bash

# Generate visualization-ready results for Patient-6 multisample (Dx + Rx)
# and hook them into the PhyloWGS Witness UI.
#
# Usage:
#   1) Make it executable:
#        chmod +x Shlush_multi-sample-test/visualize_Patient6_multisample.sh
#   2) Run:
#        ./Shlush_multi-sample-test/visualize_Patient6_multisample.sh
#   3) Start a web server from the project root, e.g.:
#        python3 -m http.server 8000
#      and open:
#        http://localhost:8000/witness/index.html
#   4) In Witness UI, select run "shlush_multi" and dataset "Patient-6_multisample"

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RESULT_DIR="${PROJECT_ROOT}/Shlush_multi-sample-test/Patient-6_multisample_chains"
WITNESS_DATA_RUN="shlush_multi"
DATASET_NAME="Patient-6_multisample"
WITNESS_DATA_DIR="${PROJECT_ROOT}/witness/data/${WITNESS_DATA_RUN}"

PYTHON2_BIN="${PYTHON2_BIN:-python2}"

echo "Project root     : ${PROJECT_ROOT}"
echo "Result directory : ${RESULT_DIR}"
echo "Witness data dir : ${WITNESS_DATA_DIR}"
echo "Dataset name     : ${DATASET_NAME}"
echo

if ! command -v "${PYTHON2_BIN}" >/dev/null 2>&1; then
  echo "ERROR: python2 not found in PATH (or via \$PYTHON2_BIN)." >&2
  echo "       Please install Python 2 or point PYTHON2_BIN to a Python 2 executable." >&2
  exit 1
fi

TREES_ZIP="${RESULT_DIR}/trees.zip"
if [[ ! -f "${TREES_ZIP}" ]]; then
  echo "ERROR: Cannot find trees.zip at ${TREES_ZIP}." >&2
  echo "       Make sure PhyloWGS multisample has completed for Patient-6." >&2
  exit 1
fi

mkdir -p "${WITNESS_DATA_DIR}"

TREES_SUMM_OUT="${WITNESS_DATA_DIR}/${DATASET_NAME}.summ.json.gz"
MUTLIST_OUT="${WITNESS_DATA_DIR}/${DATASET_NAME}.muts.json.gz"
MUTASS_OUT="${WITNESS_DATA_DIR}/${DATASET_NAME}.mutass.zip"

echo "Generating JSON / mutass results for ${DATASET_NAME}..."
echo "  Summary : ${TREES_SUMM_OUT}"
echo "  Mutlist : ${MUTLIST_OUT}"
echo "  Mutass  : ${MUTASS_OUT}"
echo

(
  cd "${PROJECT_ROOT}"
  "${PYTHON2_BIN}" write_results.py \
    "${DATASET_NAME}" \
    "${TREES_ZIP}" \
    "${TREES_SUMM_OUT}" \
    "${MUTLIST_OUT}" \
    "${MUTASS_OUT}"
)

echo
echo "Decompressing for Witness UI (browser needs plain JSON)..."
SUMM_JSON="${WITNESS_DATA_DIR}/${DATASET_NAME}.summ.json"
MUTS_JSON="${WITNESS_DATA_DIR}/${DATASET_NAME}.muts.json"
MUTASS_DIR="${WITNESS_DATA_DIR}/${DATASET_NAME}.mutass"

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
echo "  3) In the Witness UI, look for run \"${WITNESS_DATA_RUN}\" and dataset \"${DATASET_NAME}\""
echo "     to inspect trees, clusters and VAF plots."
