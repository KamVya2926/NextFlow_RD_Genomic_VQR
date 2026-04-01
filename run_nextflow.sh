#!/bin/bash
set -euo pipefail

# -----------------------------
# DNAnexus inputs (environment variables)
# -----------------------------
# DNAnexus automatically stages inputs and sets paths

NEXTFLOW_CMD="nextflow run main.nf -profile docker"

# Required inputs
NEXTFLOW_CMD+=" --samplesheet $SAMPLESHEET"
NEXTFLOW_CMD+=" --genome_file $GENOME_FILE"

# Optional file arrays
if [[ -n "${BWA_INDEX_FILES+x}" ]]; then
    for f in "${BWA_INDEX_FILES[@]}"; do
        NEXTFLOW_CMD+=" --bwa_index_files $f"
    done
fi

if [[ -n "${MINIMAP2_INDEX_FILE+x}" ]]; then
    NEXTFLOW_CMD+=" --minimap2_index_file $MINIMAP2_INDEX_FILE"
fi

if [[ -n "${GATK_REFERENCE_FILES+x}" ]]; then
    for f in "${GATK_REFERENCE_FILES[@]}"; do
        NEXTFLOW_CMD+=" --gatk_reference_files $f"
    done
fi

if [[ -n "${QSRVCFS+x}" ]]; then
    for f in "${QSRVCFS[@]}"; do
        NEXTFLOW_CMD+=" --qsrVcfs $f"
    done
fi

# Boolean options
[[ "${INDEX_GENOME:-false}" == "true" ]] && NEXTFLOW_CMD+=" --index_genome"
[[ "${FASTQC:-false}" == "true" ]] && NEXTFLOW_CMD+=" --fastqc"
[[ "${FASTP:-false}" == "true" ]] && NEXTFLOW_CMD+=" --fastp"
[[ "${BQS:-false}" == "true" ]] && NEXTFLOW_CMD+=" --bqsr"
[[ "${VARIANT_RECALIBRATION:-false}" == "true" ]] && NEXTFLOW_CMD+=" --variant_recalibration"
[[ "${DEGRADED_DNA:-false}" == "true" ]] && NEXTFLOW_CMD+=" --degraded_dna"
[[ "${IDENTITY_ANALYSIS:-false}" == "true" ]] && NEXTFLOW_CMD+=" --identity_analysis"

# Other string parameters
[[ -n "${ALIGNER:-}" ]] && NEXTFLOW_CMD+=" --aligner $ALIGNER"
[[ -n "${VARIANT_CALLER:-}" ]] && NEXTFLOW_CMD+=" --variant_caller $VARIANT_CALLER"
[[ -n "${OUTDIR:-}" ]] && NEXTFLOW_CMD+=" --outdir $OUTDIR"

# -----------------------------
# Run Nextflow
# -----------------------------
echo "Running: $NEXTFLOW_CMD"
eval $NEXTFLOW_CMD
