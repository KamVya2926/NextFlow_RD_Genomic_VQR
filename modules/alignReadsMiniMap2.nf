process minimap2Align {

    label params.platform == 'cloud' ? 'process_high' : 'process_medium'

    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-2'

    publishDir("$params.outdir/minimap2", mode: 'copy')

    input:
    tuple val(sample_id), path(reads)
    tuple path(genomeFasta), path(minimap2_index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    def r1 = reads[0]
    def r2 = reads.size() > 1 ? reads[1] : ''
    def input_reads = r2 ? "$r1 $r2" : "$r1"
    """
    echo "Running minimap2 for sample ${sample_id}"

    # Align, add read groups, and sort in one pipe
    minimap2 -t ${task.cpus} -ax sr ${minimap2_index} ${input_reads} | \
        samtools addreplacerg \
            -r "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" - | \
        samtools sort -@ 2 -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}
