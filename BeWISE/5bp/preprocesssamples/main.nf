#!/usr/bin/env nextflow

process PREPROCESS_5BP {
    publishDir params.outdir, mode: 'copy'
    
    input:
        path cx

    output:
        path "${cx}_mvalues.csv", emit: csv


    script:
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np
        

        cx_file = pd.read_csv("${cx}",header=None)
        cx_file = cx_file[cx_file.iloc[:,5] == "CG"]
        m_file = pd.DataFrame()
        m_file["chromosome"] = cx_file.iloc[:,0]
        m_file["start"] = cx_file.iloc[:,1]
        m_file["stop"] = cx_file.iloc[:,1] + 1
        sample_name = "${cx}".split("/")[-1].split(".")[0]
        m_file[sample_name] = np.log2((cx_file[3] + 100) / (cx_file[4] + 100))
        m_file.columns = ["chromosome","start","stop",sample_name]
        m_file.to_csv("${cx}_mvalues.csv",index=False)
        """
}
