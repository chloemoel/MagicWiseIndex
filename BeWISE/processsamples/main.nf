#!/usr/bin/env nextflow

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

