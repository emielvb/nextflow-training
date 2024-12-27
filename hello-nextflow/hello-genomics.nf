#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input (file of input files, one per line)
params.reads_bam = "./data/sample_bams.txt"
params.outdir = "results_genomics"
// Accessory files
params.reference        = "./data/ref/ref.fasta"
params.reference_index  = "./data/ref/ref.fasta.fai"
params.reference_dict   = "./data/ref/ref.dict"
params.intervals        = "./data/ref/intervals.bed"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink' // should not use symlink in final workflow since data would be removed on cleaning up work/ dir

    input:
        path input_bam

    output:
        tuple path(input_bam), path("${input_bam}.bai") // pair the input and output to one single channel containing both files.

    script:
    """
    samtools index '$input_bam'
    """

}

process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        tuple path(input_bam), path(input_bam_index) // paired input that corresponds to the output of SAMTOOLS_INDEX.
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.vcf"     , emit: vcf // emit uniquely names each of the outputs.
        path "${input_bam}.vcf.idx" , emit: idx
    
    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
    // note that the input_bam and input_bam_index are not given to the script, yet are still passed onto the process.
    // this is to ensure that nextflow does stage these files, such that HaplotypeCaller 
    // (which takes in these files based on naming convention, without needing them as specific arguments)
    // can access them inside the container.
}

workflow {

    // Create input channel
    reads_ch = Channel.fromPath(params.reads_bam).splitText()

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // load file paths for the accessoiry files as variables using the file() function;
    // this allows these paths to be called at will in any step in the workflow.
    ref_file = file(params.reference)
    ref_index_file = file(params.reference_index)
    ref_dict_file = file(params.reference_dict)
    intervals_file = file(params.intervals)

    // temporary diagnostics
    reads_ch.view() // the .view() function just prints the filenames to the console while running the pipeline.
    SAMTOOLS_INDEX.out.view()

    // call haplotypecaller process
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file,
    ) // note that in nextflow, all inputs to a process are positional.

}
