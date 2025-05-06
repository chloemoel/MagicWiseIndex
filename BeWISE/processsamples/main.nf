#!/usr/bin/env nextflow

process PROCESS_METHYLATION {
    publishDir params.outdir, mode: 'copy'
    
    input:
        path additional_data
        val batch_correction 
        path sample_sheet
        path sample_m_vals
    
    output:
        path "m_values_processed.csv", emit: csv

    script:
        def batch = batch_correction.collect { batch -> "\"${batch}\"" }.join(", ")
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np
        from inmoose.pycombat import pycombat_norm

        # Check if sample sheet is raw or cleaned
        with open("${sample_sheet}", "r") as s:
            first_line = s.readline().split(",")
            if first_line == [' ']:  # When it's in the raw form
                sample_info = pd.read_csv("${sample_sheet}", skiprows=8)
            else:  # Already cleaned
                sample_info = pd.read_csv("${sample_sheet}")

        # Rename sample_info columns
        sample_info.columns = ["Study_ID","Sample_Well", "Sample_Plate", "Sample_Group", "Pool_ID", "Sentrix_ID", "Sentrix_Position"]
        sample_info.dropna(subset = ["Study_ID","Sentrix_ID"], inplace=True)

        # Load m-values, samples as rows, probes as columns
        m_values = pd.read_csv("${sample_m_vals}", index_col=0)

        # Filter out probes with more than 10% missing data and fill NaN with mean of each probe
        thresh = int(len(m_values.columns) * 0.1)
        m_values.dropna(thresh=thresh, axis=0, inplace=True)
        m_values.fillna(m_values.mean(), axis=0, inplace=True)

        if "${batch_correction}" != "null":
            if "${additional_data}" != "null":
                # Load additional data and merge with sample_info
                additional_data = pd.read_csv("${additional_data}", index_col=0, dtype={0:str})
                sample_info = sample_info.merge(additional_data, left_on="Study_ID", right_index=True)

            # Create array_id for merging
            sample_info["array_id"] = sample_info["Sentrix_ID"].astype("str") + "_" + sample_info["Sentrix_Position"].astype("str")
            sample_info["array_id"] = sample_info["array_id"].apply(lambda x: x.rstrip())
            sample_info.set_index("array_id", inplace=True)

            # Merge sample info and m_values so m values are in right order 
            m_and_info = m_values.T.merge(sample_info, left_index=True, right_index=True)
            header = m_and_info["Study_ID"]
            m_values = m_and_info.drop(columns=sample_info.columns).T

            # Perform batch correction for each batch in the list
            for b in [${batch}]:
                print(m_and_info[b])
                m_values = pycombat_norm(m_values, m_and_info[b], na_cov_action="remove")

            #replace array id with study id
            m_values.columns = header

        else: #when there is no batch correction
            m_and_info = m_values.T.merge(sample_info, left_index=True, right_index=True)
            header = m_and_info["Study_ID"] 
            m_values = m_and_info.drop(columns=sample_info.columns).T
            m_values.columns = header
 
        # Save the processed m_values with samples as columns and probes as rows
        m_values.to_csv("m_values_processed.csv")
        """
}
