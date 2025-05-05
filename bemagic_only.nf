#!/usr/bin/env nextflow

include { ANNOTATE_VCF }        from './BeMAGIC/annotatevcf/main.nf'
include { CALCULATE_BEMAGIC}    from './BeMAGIC/calculatescore/main.nf'

/*
 * Pipeline parameters
 */

// Accessory files and default values
params.probe_info          = "${projectDir}/bin/data/probe_info.csv"
params.dbnsfp              = "${projectDir}/bin/data/"

workflow {

    // Create input channels from user input for sample sheet and sample directory
    vcf = Channel.fromPath(params.vcf)

    // Load the file paths for the accessory files (dbnsfp and methylation probe info)
    probe_info = file(params.probe_info)
    dbnsfp = file(params.dbnsfp)

    // Annotate vcf with VEP
    ANNOTATE_VCF(
        dbnsfp,
        vcf
    )

    // calculate BeMAGIC score
    CALCULATE_BEMAGIC(
        ANNOTATE_VCF.out,
        probe_info
    )

}
