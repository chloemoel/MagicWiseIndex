#!/usr/bin/env nextflow

include { PROCESS_METHYLATION } from './BeWISE/processsamples/main.nf'
include { CALCULATE_BEWISE }    from './BeWISE/calculatescore/main.nf'
include { CALCULATE_BEMAGIC}    from './BeMAGIC/calculatescore/main.nf'
include { ANNOTATE_VCF }        from './BeMAGIC/annotatevcf/main.nf'
include { CALCULATE_MAGICWISE } from './MagicWise/main.nf'

/*
 * Pipeline parameters
 */

// Accessory files and default values
params.probe_info          = "${projectDir}/bin/data/probe_info.csv"
params.dbnsfp              = "${projectDir}/bin/data/dbNSFP4.9a.MagicWise.txt"
params.batch_correction    = "null"
params.additional_data     = "null"
params.ens_to_gene         = "${projectDir}/bin/data/ensembl_togenename_genes_hg19.csv"

workflow {

    // Create input channels from user input for sample sheet and sample directory
    additional_data = Channel.fromPath(params.additional_data)
    batch_correction = Channel.of(params.batch_correction)
    sample_sheet = Channel.fromPath(params.sample_sheet)
    sample_m_vals= Channel.fromPath(params.sample_m_vals)
    vcf = Channel.fromPath(params.vcf)

    // Load the file paths for the accessory files (dbnsfp and methylation probe info)
    probe_info = file(params.probe_info)
    dbnsfp = file(params.dbnsfp)
    ens_to_gene = file(params.ens_to_gene)

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
    annotated_vcfs = ANNOTATE_VCF(
        dbnsfp,
        vcf
    ).collect().view()

    // calculate BeMAGIC score
    CALCULATE_BEMAGIC(
        annotated_vcfs,
        probe_info
    )

    // combine the two scores into one
    CALCULATE_MAGICWISE(
        CALCULATE_BEMAGIC.out,
        CALCULATE_BEWISE.out,
	ens_to_gene
    )
}
