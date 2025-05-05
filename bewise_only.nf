#!/usr/bin/env nextflow

include { PROCESS_METHYLATION } from './BeWISE/processsamples/main.nf'
include { CALCULATE_BEWISE }    from './BeWISE/calculatescore/main.nf'

/*
 * Pipeline parameters
 */

// Accessory files and default values
params.probe_info          = "${projectDir}/bin/data/probe_info.csv"
params.batch_correction    = "null"
params.additional_data     = "null"

workflow {

    // Create input channels from user input for sample sheet and sample directory
    additional_data = Channel.fromPath(params.additional_data)
    batch_correction = Channel.of(params.batch_correction)
    sample_sheet = Channel.fromPath(params.sample_sheet)
    sample_m_vals= Channel.fromPath(params.sample_m_vals)

    // Load the file paths for the accessory files (dbnsfp and methylation probe info)
    probe_info = file(params.probe_info)

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
}
