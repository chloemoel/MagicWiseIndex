#!/usr/bin/env nextflow

include { PROCESS_METHYLATION }         from './BeWISE/processsamples/main.nf'
include { CALCULATE_BEWISE }         from './BeWISE/calculatescore/main.nf'
include { CALCULATE_BEMAGIC}    from './BeMAGIC/calculatescore/main.nf'
include { ANNOTATE_VCF }            from './BeMAGIC/annotatevcf/main.nf'
include { CALCULATE_MAGICWISE }                 from './MagicWise/main.nf'

/*
 * Pipeline parameters
 */

// Accessory files and default values
params.probe_info          = "${projectDir}/data/probe_info.csv"
params.cache_dir           = "${projectDir}/data/vep_data"
params.genome_fasta        = "${projectDir}/data/vep_data/homo_sapiens/113_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
params.cadd                = "${projectDir}/data/vep_data/CADD/whole_genome_SNVs.tsv.gz"
params.cadd_tbi            = "${projectDir}/data/vep_data/CADD/whole_genome_SNVs.tsv.gz.tbi"
params.batch_correction    = "null"
params.additional_data     = "null"
params.vcf                 = "null"

workflow {

    // Create input channels from user input for sample sheet and sample directory
    additional_data = Channel.fromPath(params.additional_data)
    batch_correction = Channel.of(params.batch_correction)
    sample_sheet = Channel.fromPath(params.sample_sheet)
    sample_m_vals= Channel.fromPath(params.sample_m_vals)
    vcf = Channel.fromPath(params.vcf)

    // Load the file paths for the accessory files (reference and intervals)
    probe_info = file(params.probe_info)
    cache_dir = file(params.cache_dir)
    genome_fasta = file(params.genome_fasta)
    cadd = file(params.cadd)
    cadd_tbi = file(params.cadd_tbi)

    //Clean up data and assess for batch correction
    PROCESS_METHYLATION(
        additional_data,
        batch_correction,
        sample_sheet,
        sample_m_vals
    )

    // Calculate BeWISE score
    CALCULATE_BEWISE(
        PROCESS_METHYLATION.out,
        probe_info
    )

    // Annotate vcf with VEP
    ANNOTATE_VCF(
        cache_dir,
        vcf,
        genome_fasta,
        cadd,
        cadd_tbi
    )

    // calculate BeMAGIC score
    CALCULATE_BEMAGIC(
        ANNOTATE_VCF.out,
        probe_info
    )

    // combine the two scores into one
    CALCULATE_MAGICWISE(
        CALCULATE_BEMAGIC.out,
        CALCULATE_BEWISE.out
    )
}
