// Use newest nextflow dsl
nextflow.enable.dsl = 2

log.info """\
============================================
      DNASeq Pipeline Configuration
============================================
platform             : ${params.platform}
samplesheet          : ${params.samplesheet}
genome               : ${params.genome_file}
bwa index            : ${params.bwa_index_files}
minimap2 index       : ${params.minimap2_index_file}
gatk reference       : ${params.gatk_reference_files}
index genome         : ${params.index_genome}
qsr truth vcfs       : ${params.qsrVcfs}
output directory     : ${params.outdir}
fastqc               : ${params.fastqc}
aligner              : ${params.aligner}
variant caller       : ${params.variant_caller}
bqsr                 : ${params.bqsr}
degraded_dna         : ${params.degraded_dna}
variant_recalibration: ${params.variant_recalibration}
identity_analysis    : ${params.identity_analysis}
============================================
""".stripIndent()

// -----------------------------
// Module includes (unchanged)
// -----------------------------
if (params.index_genome) include { indexGenome } from './modules/indexGenome'
if (params.fastqc) include { FASTQC } from './modules/FASTQC'
if (params.fastp) include { fastpReads } from './modules/fastp'
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
if (params.bqsr) include { baseRecalibrator } from './modules/BQSR'
include { combineGVCFs; genotypeGVCFs } from './modules/processGVCFs'
if (params.variant_recalibration) {
    include { variantRecalibrator } from './modules/variantRecalibrator'
} else {
    include { filterVCF } from './modules/filterVCF'
}
if (params.identity_analysis) include { identityAnalysis } from './modules/identityAnalysis'

if (params.aligner == 'bwa-mem') {
    include { alignReadsBwaMem } from './modules/alignReadsBwaMem'
} else if (params.aligner == 'bwa-aln') {
    include { alignReadsBwaAln } from './modules/alignReadsBwaAln'
} else if (params.aligner == 'minimap2') {
    include { minimap2Align } from './modules/alignReadsMiniMap2'
    include { indexGenomeMinimap2 } from './modules/indexGenomeMinimap2'
} else {
    error "Unsupported aligner: ${params.aligner}"
}

if (params.variant_caller == 'haplotype-caller') {
    include { haplotypeCaller } from './modules/haplotypeCaller'
} else {
    error "Unsupported variant caller"
}

if (params.degraded_dna) {
    include { mapDamage2 } from './modules/mapDamage'
    include { indexMapDamageBam } from './modules/indexBam'
}

// -----------------------------
// WORKFLOW
// -----------------------------
workflow {

    // ✅ FIX: arrays → Channel.from()
    qsrc_vcf_ch = Channel.from(params.qsrVcfs)

    // samplesheet stays fromPath ✅
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row ->
            row.size() == 4 ? tuple(row[0], [row[1], row[2]]) :
            row.size() == 3 ? tuple(row[0], [row[1]]) :
            error "Bad samplesheet row: $row"
        }

    if (params.fastqc) FASTQC(read_pairs_ch)

    preprocessed_reads_ch = params.fastp ?
        fastpReads(read_pairs_ch).map { id, r1, r2, _, _ -> tuple(id, [r1, r2]) } :
        read_pairs_ch

    // -----------------------------
    // ALIGNMENT
    // -----------------------------
    if (params.aligner in ['bwa-mem','bwa-aln']) {

        reference_genome_ch = params.index_genome ?
            indexGenome(params.genome_file).flatten().collect() :
            Channel.from(params.bwa_index_files).collect()   // ✅ FIX

        align_ch = (params.aligner == 'bwa-mem') ?
            alignReadsBwaMem(preprocessed_reads_ch, reference_genome_ch) :
            alignReadsBwaAln(preprocessed_reads_ch, reference_genome_ch)

    } else if (params.aligner == 'minimap2') {

        minimap2_index_ch = params.index_genome ?
            indexGenomeMinimap2(params.genome_file) :
            Channel.of(tuple(params.genome_file, params.minimap2_index_file)) // ✅ FIX (removed file())

        align_ch = minimap2Align(preprocessed_reads_ch, minimap2_index_ch)

        reference_genome_ch = params.index_genome ?
            indexGenome(params.genome_file).flatten().collect() :
            Channel.from(params.gatk_reference_files).collect() // ✅ FIX
    }

    // -----------------------------
    // POST-ALIGNMENT
    // -----------------------------
    sort_ch = sortBam(align_ch)
    mark_ch = markDuplicates(sort_ch)
    indexed_bam_ch = indexBam(mark_ch)

    mapDamage_ch = params.degraded_dna ?
        indexMapDamageBam(mapDamage2(indexed_bam_ch, reference_genome_ch)) :
        indexed_bam_ch

    // ✅ FIX: arrays → Channel.from()
    knownSites_ch = Channel.from(params.qsrVcfs)
        .filter { it.name.endsWith('.tbi') || it.name.endsWith('.idx') }
        .map { "--known-sites ${it}" }
        .collect()

    bqsr_ch = params.bqsr ?
        baseRecalibrator(mapDamage_ch, knownSites_ch, reference_genome_ch, qsrc_vcf_ch.collect()) :
        mapDamage_ch

    // -----------------------------
    // VARIANT CALLING
    // -----------------------------
    gvcf_ch = haplotypeCaller(bqsr_ch, reference_genome_ch).collect()

    all_gvcf_ch = gvcf_ch.collect { list ->
        tuple(
            list.collate(3)*.getAt(0),
            list.collate(3)*.getAt(1),
            list.collate(3)*.getAt(2)
        )
    }

    combined_gvcf_ch = combineGVCFs(all_gvcf_ch, reference_genome_ch)
    final_vcf_ch     = genotypeGVCFs(combined_gvcf_ch, reference_genome_ch)

    // -----------------------------
    // FILTERING
    // -----------------------------
    if (params.variant_recalibration) {

        knownSitesArgs_ch = Channel.from(params.qsrVcfs) // ✅ FIX
            .filter { it.name.endsWith('.vcf') || it.name.endsWith('.vcf.gz') }
            .map { "--resource:${it.baseName} ${it}" }
            .collect()

        filtered_vcf_ch = variantRecalibrator(
            final_vcf_ch,
            knownSitesArgs_ch,
            reference_genome_ch,
            qsrc_vcf_ch.collect()
        )

    } else {
        filtered_vcf_ch = filterVCF(final_vcf_ch, reference_genome_ch)
    }

    // -----------------------------
    // IDENTITY ANALYSIS (FIXED TMP)
    // -----------------------------
    if (params.identity_analysis) {

        psam_info_ch = Channel.fromPath(params.samplesheet)
            .splitCsv(sep: '\t')
            .map { row ->
                tuple(row[0], row.size()==4 ? row[3] : row[2])
            }

        psam_file_ch = psam_info_ch.map { sample_id, sex ->
            def f = new File("${workflow.workDir}/combined_samples.psam") // ✅ FIX
            f.text = "#IID\tSID\tPAT\tMAT\tSEX\n${sample_id}\t${sample_id}\t0\t0\t${sex ?: 'NA'}\n"
            f
        }

        identityAnalysis(filtered_vcf_ch, psam_file_ch)
    }
}

workflow.onComplete {
    log.info(workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong")
}
