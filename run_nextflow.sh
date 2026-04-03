#!/bin/bash
set -euo pipefail

echo "=== DEBUG: Listing inputs ==="
env | sort | grep -E 'sample|genome|bwa|minimap|gatk|qsr' || true
echo "============================"

# -----------------------------
# Validate required inputs
# -----------------------------
if [[ -z "${samplesheet:-}" ]]; then
    echo "ERROR: samplesheet is not set"
    exit 1
fi

if [[ -z "${genome_file:-}" ]]; then
    echo "ERROR: genome_file is not set"
    exit 1
fi

echo "Samplesheet: $samplesheet"
echo "Genome: $genome_file"

echo "=== Files in working dir ==="
ls -lh
echo "============================"

# -----------------------------
# Build Nextflow command SAFELY
# -----------------------------
cmd=(nextflow run main.nf -profile docker)

# Required
cmd+=(--samplesheet "$samplesheet")
cmd+=(--genome_file "$genome_file")

# -----------------------------
# Optional arrays
# -----------------------------

# IMPORTANT: DNAnexus may pass arrays as space-separated strings
# So we normalize them

if [[ -n "${bwa_index_files:-}" ]]; then
    echo "Using BWA index files: $bwa_index_files"
    for f in $bwa_index_files; do
        cmd+=(--bwa_index_files "$f")
    done
fi

if [[ -n "${gatk_reference_files:-}" ]]; then
    echo "Using GATK reference files: $gatk_reference_files"
    for f in $gatk_reference_files; do
        cmd+=(--gatk_reference_files "$f")
    done
fi

if [[ -n "${qsrVcfs:-}" ]]; then
    echo "Using QSR VCFs: $qsrVcfs"
    for f in $qsrVcfs; do
        cmd+=(--qsrVcfs "$f")
    done
fi

# -----------------------------
# Optional single files
# -----------------------------
if [[ -n "${minimap2_index_file:-}" ]]; then
    cmd+=(--minimap2_index_file "$minimap2_index_file")
fi

# -----------------------------
# Booleans
# -----------------------------
[[ "${index_genome:-false}" == "true" ]] && cmd+=(--index_genome)
[[ "${fastqc:-false}" == "true" ]] && cmd+=(--fastqc)
[[ "${fastp:-false}" == "true" ]] && cmd+=(--fastp)
[[ "${bqsr:-false}" == "true" ]] && cmd+=(--bqsr)
[[ "${variant_recalibration:-false}" == "true" ]] && cmd+=(--variant_recalibration)
[[ "${degraded_dna:-false}" == "true" ]] && cmd+=(--degraded_dna)
[[ "${identity_analysis:-false}" == "true" ]] && cmd+=(--identity_analysis)

# -----------------------------
# Strings
# -----------------------------
[[ -n "${aligner:-}" ]] && cmd+=(--aligner "$aligner")
[[ -n "${variant_caller:-}" ]] && cmd+=(--variant_caller "$variant_caller")
[[ -n "${outdir:-}" ]] && cmd+=(--outdir "$outdir")

# -----------------------------
# Run
# -----------------------------
echo "=== Running Nextflow ==="
printf '%q ' "${cmd[@]}"
echo
echo "========================"

"${cmd[@]}"
