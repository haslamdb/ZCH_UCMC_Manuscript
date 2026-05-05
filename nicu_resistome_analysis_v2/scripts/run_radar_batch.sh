#!/bin/bash
# Parallel RADAR batch via biogpu unified_pipeline.
# Splits the sample manifest into N chunks (one per GPU) and runs unified_pipeline
# concurrently across all available GPUs. Each instance processes its chunk
# sample-by-sample but uses a single orchestrator + DIAMOND + bloom prescreen,
# so the DB load amortization comes from staying inside one orchestrator process
# and from DIAMOND's startup being trivial.

set -euo pipefail

REPO_ROOT="/home/david/projects/biogpu"
V2="/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis_v2"
MANIFEST="${V2}/data/sample_manifest.csv"
OUT_DIR="${V2}/radar_runs/unified"
LOG_DIR="${V2}/logs/batch"
CHUNK_DIR="${LOG_DIR}/chunks"
mkdir -p "${OUT_DIR}" "${LOG_DIR}" "${CHUNK_DIR}"

UNI="${REPO_ROOT}/build/runtime/kernels/unified_pipeline"
REF_DB="${REPO_ROOT}/data/fq_resistance_db_v2/nucleotide"
RES_DB="${REPO_ROOT}/data/fq_resistance_db_v2/protein"
BLOOM="${REPO_ROOT}/data/amr_protein_db/protein_bloom.bin"
DIAMOND="${REPO_ROOT}/data/amr_protein_db/proteins"
META="${REPO_ROOT}/data/amr_protein_db/protein_details.json"

# GPUs to use (PCI ordering: 0=Blackwell, 1=A6000, 2=Blackwell)
GPUS=(0 1 2)
NGPU=${#GPUS[@]}

# Build per-GPU CSVs in unified_pipeline format: SampleName,FilePath,R1 file,R2 file
HEADER="SampleName,FilePath,R1 file,R2 file"
for i in "${!GPUS[@]}"; do
  echo "${HEADER}" > "${CHUNK_DIR}/chunk_${GPUS[$i]}.csv"
done

# Skip already-completed samples (gene_stats.tsv + allele_frequencies.tsv both present)
TOTAL=0
ASSIGNED=0
SKIPPED=0
i=0
while IFS=, read -r sid r1 r2; do
  [[ -z "$sid" || "$sid" == "sample_id" ]] && continue
  TOTAL=$((TOTAL + 1))
  if [[ -f "${OUT_DIR}/${sid}/${sid}_gene_stats.tsv" ]] \
     && [[ -f "${OUT_DIR}/${sid}/${sid}_allele_frequencies.tsv" ]]; then
    SKIPPED=$((SKIPPED + 1))
    continue
  fi
  gpu_idx=$((ASSIGNED % NGPU))
  gpu="${GPUS[$gpu_idx]}"
  fpath=$(dirname "$r1")
  r1f=$(basename "$r1")
  r2f=$(basename "$r2")
  echo "${sid},${fpath},${r1f},${r2f}" >> "${CHUNK_DIR}/chunk_${gpu}.csv"
  ASSIGNED=$((ASSIGNED + 1))
done < "${MANIFEST}"

echo "manifest=${TOTAL} skipped=${SKIPPED} assigned=${ASSIGNED}"
for g in "${GPUS[@]}"; do
  n=$(($(wc -l < "${CHUNK_DIR}/chunk_${g}.csv") - 1))
  echo "  gpu${g}: ${n} samples"
done
[[ $ASSIGNED -eq 0 ]] && { echo "Nothing to do."; exit 0; }

run_chunk() {
  local gpu=$1
  local chunk="${CHUNK_DIR}/chunk_${gpu}.csv"
  local log="${LOG_DIR}/worker_gpu${gpu}.log"
  local n=$(($(wc -l < "${chunk}") - 1))
  if [[ $n -le 0 ]]; then
    echo "[gpu${gpu}] no samples"
    return 0
  fi
  echo "[gpu${gpu}] starting ${n} samples"
  cd "${REPO_ROOT}"
  CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=${gpu} "${UNI}" \
    --csv "${chunk}" \
    --output-dir "${OUT_DIR}" \
    --reference-db "${REF_DB}" \
    --resistance-db "${RES_DB}" \
    --protein-bloom "${BLOOM}" \
    --diamond-db "${DIAMOND}" \
    --amr-metadata "${META}" \
    --gpu-device 0 \
    > "${log}" 2>&1 || echo "[gpu${gpu}] worker exited non-zero (some samples may still have outputs)"
  echo "[gpu${gpu}] done"
}

PIDS=()
for g in "${GPUS[@]}"; do
  run_chunk "${g}" &
  PIDS+=($!)
done
echo "PIDs: ${PIDS[@]}"
echo "Follow:"
for g in "${GPUS[@]}"; do echo "  tail -f ${LOG_DIR}/worker_gpu${g}.log"; done
wait "${PIDS[@]}" || true
echo "all workers finished"

# Quick completeness count
COMPLETED=0
for d in "${OUT_DIR}"/*/; do
  sid=$(basename "$d")
  if [[ -f "${d}/${sid}_gene_stats.tsv" ]] && [[ -f "${d}/${sid}_allele_frequencies.tsv" ]]; then
    COMPLETED=$((COMPLETED + 1))
  fi
done
echo "completed: ${COMPLETED}/${TOTAL}"

# Post-hoc aggregation: parse pipe-delim protein_id and join ARO + MEGARes
# annotations from data/amr_protein_db/amr_hierarchy_annotated.tsv. Required
# because the in-pipeline coverage_calculator expects 11-part DIAMOND headers
# but the current DB uses the compact 4-part form.
AGG="${V2}/scripts/aggregate_radar_outputs.py"
if [[ -f "${AGG}" ]]; then
  echo "running aggregator: ${AGG}"
  python3 "${AGG}" --v2-root "${V2}" --biogpu-root "${REPO_ROOT}" \
    > "${LOG_DIR}/aggregate.log" 2>&1 \
    && echo "aggregation done -> ${V2}/data/" \
    || echo "[WARN] aggregator failed; see ${LOG_DIR}/aggregate.log"
else
  echo "[WARN] aggregator not found at ${AGG}; skipping"
fi
