#!/usr/bin/env nextflow

include { PROCESS_SAMPLES }         from './BeWISE/processsamples/main.nf'
include { CALCULATE_SCORE }         from './BeWISE/calculatescore/main.nf'
include { CALCULATE_GENE_SCORE }    from './Genetic/calculatescore/main.nf'
include { ANNOTATE_VCF }            from './Genetic/annotatevcf/main.nf'
include { COMBINE }                 from './Combine/main.nf'

/*
 * Pipeline parameters
 */

// Accessory files and default values
params.probe_info          = "${projectDir}/data/ref/probe_info.csv"
params.genome_fasta        = "${projectDir}/data/ref/GRCh37_genome.fa"
params.CADD                = "${projectDir}/test_data"
params.batch_correction    = null
params.additional_data     = null

workflow {

    // Create input channels from user input for sample sheet and sample directory
    additional_data = Channel.fromPath(params.additional_data)
    batch_correction = Channel.of(params.batch_correction)
    sample_sheet = Channel.fromPath(params.sample_sheet)
    sample_m_vals= Channel.fromPath(params.sample_m_vals)
    gvcf = Channel.fromPath(params.gvcf)
    CADD = Channel.fromPath(params.CADD)
    

    // Load the file paths for the accessory files (reference and intervals)
    genome_fasta = file(params.genome_fasta)
    probe_info = file(params.probe_info)

    //Clean up data and assess for batch correction
    PROCESS_SAMPLES(
        additional_data,
        batch_correction,
        sample_sheet,
        sample_m_vals
    )

    // Calculate BeWISE score
    CALCULATE_SCORE(
        PROCESS_SAMPLES.out,
        probe_info
    )

    // Annotate gvcf with VEP
    ANNOTATE_VCF(
        genome_fasta,
        gvcf,
        CADD
    )

    // calculate genetic score
    CALCULATE_GENE_SCORE(
        ANNOTATE_VCF.out,
        probe_info,
        sample_sheet
    )

    // combine the two scores into one
    COMBINE(
        CALCULATE_GENE_SCORE.out,
        CALCULATE_SCORE.out
    )
}
