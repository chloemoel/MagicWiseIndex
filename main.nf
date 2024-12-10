#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Accessory files and default values
params.probe_info          = "${projectDir}/data/ref/probe_info.csv"
params.batch_correction    = null
params.additional_data     = null

process PROCESS_SAMPLES {

    conda /*TODO add conda info*/

    publishDir 'BeWISE_outputs', mode: 'symlink'

    input:
        path additional_data
        val batch_correction
        path sample_sheet
        path sample_m_vals
        
    output:
        path("${sample_m_vals}_processed.csv"), emit: csv

    script:
    '''
    python3 data_cleanup.py \
        -a ${additional_data} \
        -b ${batch_correction} \
        ${sample_sheet} \
        ${sample_m_vals}
    '''
}

process CALCULATE_SCORE {

    conda /*ADD CONDA*/

    publishDir 'BeWISE_outputs', mode: 'symlink'

    input:
       path sample_m_values 
       path probe_info
    
    output:
        path "bewise_scores.csv", emit: csv

    script:
    '''
    python3 calculate_score.py \
        ${probe_info} \
        ${sample_m_vals}
    '''

}


workflow {

    // Create input channels from user input for sample sheet and sample directory
    additional_data = Channel.fromPath(params.additional_data)
    batch_correction = Channel.of(params.batch_correction)
    sample_sheet = Channel.fromPath(params.sample_sheet)
    sample_m_vals= Channel.fromPath(params.sample_m_vals)
    

    // Load the file paths for the accessory files (reference and intervals)
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
}
