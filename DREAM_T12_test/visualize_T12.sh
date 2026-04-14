#!/usr/bin/env bash

# Prepare PhyloWGS Witness visualization for both DREAM T12 runs:
#   - T12_raw            (from raw/)
#   - T12_with_CN_limit  (from with_CN_limit/)
#
# Usage:
#   1) chmod +x DREAM_T12_test/visualize_T12.sh
#   2) ./DREAM_T12_test/visualize_T12.sh
#   3) From the project root, start a web server:
#        python3 -m http.server 8000
#   4) Open: http://localhost:8000/witness/index.html
#   5) In the Witness UI, select run "dream_T12" to see both datasets.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PYTHON2_BIN="${PYTHON2_BIN:-python2}"

WITNESS_RUN="dream_T12"
WITNESS_DATA_DIR="${PROJECT_ROOT}/witness/data/${WITNESS_RUN}"

if ! command -v "${PYTHON2_BIN}" >/dev/null 2>&1; then
    echo "ERROR: python2 not found. Set \$PYTHON2_BIN if needed." >&2
    exit 1
fi

prepare_dataset() {
    local result_dir="$1"
    local dataset_name="$2"
    local trees_zip="${result_dir}/chains/trees.zip"

    echo "============================================"
    echo "  Preparing: ${dataset_name}"
    echo "  Source:    ${result_dir}"
    echo "============================================"

    if [[ ! -f "${trees_zip}" ]]; then
        echo "ERROR: ${trees_zip} not found. Run run_phylowgs.sh first." >&2
        return 1
    fi

    mkdir -p "${WITNESS_DATA_DIR}"

    local summ_gz="${WITNESS_DATA_DIR}/${dataset_name}.summ.json.gz"
    local muts_gz="${WITNESS_DATA_DIR}/${dataset_name}.muts.json.gz"
    local mutass_zip="${WITNESS_DATA_DIR}/${dataset_name}.mutass.zip"

    # If write_results outputs already exist in result_dir, reuse them;
    # otherwise generate from trees.zip.
    local existing_summ="${result_dir}/${dataset_name}.summ.json.gz"
    local existing_muts="${result_dir}/${dataset_name}.muts.json.gz"
    local existing_mutass="${result_dir}/${dataset_name}.mutass.zip"

    if [[ -f "${existing_summ}" && -f "${existing_muts}" && -f "${existing_mutass}" ]]; then
        echo "  Reusing existing write_results output from ${result_dir}"
        cp "${existing_summ}" "${summ_gz}"
        cp "${existing_muts}" "${muts_gz}"
        cp "${existing_mutass}" "${mutass_zip}"
    else
        echo "  Generating results with write_results.py..."
        (
            cd "${PROJECT_ROOT}"
            "${PYTHON2_BIN}" write_results.py \
                "${dataset_name}" \
                "${trees_zip}" \
                "${summ_gz}" \
                "${muts_gz}" \
                "${mutass_zip}"
        )
    fi

    echo "  Decompressing for Witness UI..."
    gunzip -c "${summ_gz}" > "${WITNESS_DATA_DIR}/${dataset_name}.summ.json"
    gunzip -c "${muts_gz}" > "${WITNESS_DATA_DIR}/${dataset_name}.muts.json"

    local mutass_dir="${WITNESS_DATA_DIR}/${dataset_name}.mutass"
    rm -rf "${mutass_dir}"
    mkdir -p "${mutass_dir}"
    unzip -o -j "${mutass_zip}" -d "${mutass_dir}" >/dev/null 2>&1

    echo "  Done: ${dataset_name}"
    echo ""
}

prepare_dataset "${SCRIPT_DIR}/raw"           "T12_raw"
prepare_dataset "${SCRIPT_DIR}/with_CN_limit" "T12_with_CN_limit"

echo "Updating Witness index..."
(
    cd "${PROJECT_ROOT}/witness"
    "${PYTHON2_BIN}" index_data.py
)

echo ""
echo "All done. Next steps:"
echo "  1) cd ${PROJECT_ROOT}"
echo "  2) python3 -m http.server 8000"
echo "  3) Open http://localhost:8000/witness/index.html"
echo "  4) Select run \"${WITNESS_RUN}\" → datasets \"T12_raw\" and \"T12_with_CN_limit\""
