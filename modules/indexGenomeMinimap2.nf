process indexGenomeMinimap2 {

    label params.platform == 'cloud' ? 'process_medium' : 'process_low'

    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:e1ea28074233d7265a5dc2111d6e55130dff5653-2'

    publishDir("$params.outdir/GENOME_IDX", mode: 'copy', saveAs: { fn -> fn.endsWith('.mmi') ? fn : null })

    input:
    path genomeFasta

    output:
    tuple path(genomeFasta), path("${genomeFasta}.mmi")

    script:
    """
    echo "Building minimap2 index for ${genomeFasta}"
    minimap2 -d ${genomeFasta}.mmi ${genomeFasta}
    echo "Minimap2 indexing done."
    """
}
