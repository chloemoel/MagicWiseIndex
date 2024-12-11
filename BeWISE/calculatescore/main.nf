#!/usr/bin/env nextflow

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

