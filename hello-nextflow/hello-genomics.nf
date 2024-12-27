#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
params.reads_bam = "./data/bam/reads_mother.bam"
params.outdir = "results_genomics"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink' // should not use symlink in final workflow since data would be removed on cleaning up work/ dir

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    script:
    """
    samtools index '$input_bam'
    """

}

workflow {

    // Create input channel
    reads_ch = Channel.fromPath(params.reads_bam)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

}
