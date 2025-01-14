#!/usr/bin/env nextflow

process PROCESS_SAMPLES {
    container "oras://community.wave.seqera.io/library/inmoose_pip_numpy_pandas:2c57c3780755c264"
    
    publishDir params.outdir, mode: 'symlink'

    input:
        path additional_data
        val batch_correction
        path sample_sheet
        path sample_m_vals
        
    output:
        path("m_values_processed.csv"), emit: csv

    script:
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np
        from inmoose.pycombat import pycombat_norm

        m_values = pd.read_csv("${sample_m_vals}",index_col=0)

        ## Check to see if the first lines are from the sample sheet or if already cleaned.
        with open("${sample_sheet}", "r") as s:
            if s.readline().split(",") == " ":
                sample_info = pd.read_csv("${sample_sheet}", skiprows=7,index_col=0)
            else:
                sample_info = pd.read_csv("${sample_sheet}", index_col=0)

        # filter out sites with more than 10% of samples missing that probe. then fill NA with probe mean
        thresh = int(len(m_values) * 0.1)
        m_values.dropna(thresh=thresh , inplace=True, axis = 1)

        #fills NA values with the mean of the probe
        m_values.fillna(m_values.mean(), axis=0, inplace=True)

        if "${batch_correction}" != "null":
            if "${additional_data}" != "null":
                additional_data = pd.read_csv("${additional_data}", index_col=0)
                sample_info = sample_info.merge(additional_data, left_index=True, right_index=True)

            sample_info["array_id"] = sample_info["Sentrix_ID"].astype("str") + "_" + sample_info["Sentrix_Position"].astype("str")

            sample_info.sort_values(by = "array_id", key = m_values.columns, axis=0, inplace=True)
            batch_columns = np.array(${batch_correction})

            #perform batch correction
            for b in batch_columns:
                m_values = pycombat_norm(m_values.drop(columns=batch_columns), m_values[b])
                print("Batch correction for", b ,"done")

        m_values.to_csv("m_values_processed.csv")
        """
}

