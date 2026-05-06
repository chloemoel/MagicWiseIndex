#!/usr/bin/env nextflow

include { PROCESS_METHYLATION_ARRAY } from './BeWISE/array/processsamples/main.nf'
include { CALCULATE_BEWISE_ARRAY }    from './BeWISE/array/calculatescore/main.nf'
include { PROCESS_METHYLATION_BIS } from './BeWISE/bisulfite/processsamples/main.nf'
include { CALCULATE_BEWISE_BIS }    from './BeWISE/bisulfite/calculatescore/main.nf'

/*
 * Pipeline parameters
 */


// Accessory files and default values
probe_info		= file("${projectDir}/bin/data/probe_info.csv")


workflow {


 if (params.sequencing_type == "bisulfite") {
    batch_correction	= params.batch_correction	? Channel.fromPath(params.batch_correction) : Channel.value([])
    additional_data		= params.additional_data	? Channel.fromPath(params.additional_data) : Channel.value([])
    bisulfite_files 	= Channel.fromPath(params.bisulfite_files).collect()


        PROCESS_METHYLATION_BIS(
            batch_correction,
            bisulfite_files,
            additional_data,
        )

        CALCULATE_BEWISE_BIS(
            PROCESS_METHYLATION_BIS.out,
            probe_info
        )

    } else {

        batch_correction	= params.batch_correction	? Channel.fromPath(params.batch_correction) : Channel.value([])
        additional_data		= params.additional_data	? Channel.fromPath(params.additional_data) : Channel.value([])
        sample_sheet		= Channel.fromPath(params.sample_sheet) 
        sample_m_vals 		= Channel.fromPath(params.sample_m_vals) 

        PROCESS_METHYLATION_ARRAY(
            batch_correction,
            sample_sheet,
            sample_m_vals,
            additional_data
        )

        CALCULATE_BEWISE_ARRAY(
            PROCESS_METHYLATION_ARRAY.out,
            probe_info
        )
    }

}
