#!/bin/bash
set -euo pipefail

PHYLOWGS_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON=python2

CNV_FILE="${DATA_DIR}/T12_CN.txt"
VCF_RAW="${DATA_DIR}/T12-noXY_subsampled.vcf"
VCF_CN_FILTER="${DATA_DIR}/T12-noXY_subsampled_with_CN_filter.vcf"
OUT_RAW="${DATA_DIR}/raw"
OUT_CN_LIMIT="${DATA_DIR}/with_CN_limit"

NUM_CHAINS=4
BURNIN=1000
MCMC_SAMPLES=2500
MH_ITERATIONS=5000

run_phylowgs() {
    local vcf_file="$1"
    local output_dir="$2"
    local run_name="$3"

    echo "============================================"
    echo "  Running PhyloWGS: ${run_name}"
    echo "  VCF:  ${vcf_file}"
    echo "  CNV:  ${CNV_FILE}"
    echo "  Output: ${output_dir}"
    echo "============================================"

    mkdir -p "${output_dir}"

    local ssm_file="${output_dir}/ssm_data.txt"
    local cnv_file="${output_dir}/cnv_data.txt"
    local params_file="${output_dir}/params.json"

    echo "[Step 1/3] Parsing inputs..."
    ${PYTHON} "${PHYLOWGS_DIR}/parser/create_phylowgs_inputs.py" \
        --vcf-type "T12=mutect_smchet" \
        --cnvs "T12=${CNV_FILE}" \
        --tumor-sample tumor \
        --output-variants "${ssm_file}" \
        --output-cnvs "${cnv_file}" \
        --output-params "${params_file}" \
        "T12=${vcf_file}"

    echo "[Step 2/3] Running MCMC (multievolve, ${NUM_CHAINS} chains)..."
    ${PYTHON} "${PHYLOWGS_DIR}/multievolve.py" \
        --num-chains ${NUM_CHAINS} \
        --ssms "${ssm_file}" \
        --cnvs "${cnv_file}" \
        -O "${output_dir}/chains" \
        -p "${params_file}" \
        -B ${BURNIN} \
        -s ${MCMC_SAMPLES} \
        -i ${MH_ITERATIONS}

    echo "[Step 3/3] Writing results..."
    ${PYTHON} "${PHYLOWGS_DIR}/write_results.py" \
        --include-ssm-names \
        "${run_name}" \
        "${output_dir}/chains/trees.zip" \
        "${output_dir}/${run_name}.summ.json.gz" \
        "${output_dir}/${run_name}.muts.json.gz" \
        "${output_dir}/${run_name}.mutass.zip"

    echo "Done: ${run_name}"
    echo ""
}

run_phylowgs "${VCF_RAW}"       "${OUT_RAW}"       "T12_raw"
run_phylowgs "${VCF_CN_FILTER}" "${OUT_CN_LIMIT}"  "T12_with_CN_limit"

echo "All runs complete."
echo "Results:"
echo "  Raw:          ${OUT_RAW}/"
echo "  With CN limit: ${OUT_CN_LIMIT}/"
