nextflow.enable.dsl=2

process fastpReads {

    tag "$sample_id"

    container 'staphb/fastp:latest'

    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
tuple val(sample_id),
      path("${sample_id}_R1.fastq.gz"),
      path("${sample_id}_R2.fastq.gz"),
      path("${sample_id}_fastp.html"),
      path("${sample_id}_fastp.json")
    when:
    params.fastp

    script:
    def is_paired = reads.size() == 2

    if (is_paired) {
        """
        fastp \
            -i ${reads[0]} -I ${reads[1]} \
            -o ${sample_id}_R1.fastq.gz \
            -O ${sample_id}_R2.fastq.gz \
            --thread 4 \
            --html ${sample_id}_fastp.html \
            --json ${sample_id}_fastp.json
        """
    }
    else {
        """
        fastp \
            -i ${reads[0]} \
            -o ${sample_id}_R1.fastq.gz \
            --thread 4 \
            --html ${sample_id}_fastp.html \
            --json ${sample_id}_fastp.json

        touch ${sample_id}_R2.fastq.gz
        """
    }
}
