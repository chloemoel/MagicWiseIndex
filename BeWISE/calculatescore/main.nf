#!/usr/bin/env nextflow

process CALCULATE_SCORE {
    container "oras://community.wave.seqera.io/library/inmoose_pip_numpy_pandas:2c57c3780755c264"

    publishDir params.outdir, mode: 'symlink'

    input:
       path processed_m_values 
       path probe_info
    
    output:
        path "BeWISE_scores.csv", emit: csv

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np

    probe_info = pd.read_csv("${probe_info}", index_col=0)
    m_values = pd.read_csv("${processed_m_values}", index_col=0)

    ## Merge m values with probe information
    score = m_values.merge(probe_info[["weights"]], left_index=True, right_index=True)

    sample_columns = score.drop(columns=["weights"])
    cpg_weights = score["weights"]

    # Multiply cpg values with cpg_weights
    weighted_cpgs = sample_columns.multiply(cpg_weights, axis="index")

    #this merges our weighted cpg site with the info on the gene
    full = weighted_cpgs.merge(probe_info, left_index=True, right_on="ID")
    full.drop(columns=["weights","Gene"], inplace=True)
    full.drop_duplicates(inplace=True)

    # This collapses the cpgs (that are weighted by chromatin state) into a single value per gene
    collapse = full.groupby(["Symbol", "gc_content","gene_length", "Transcript count"]).agg('sum')
    final = collapse.reset_index()
    final.index = final["Symbol"]
    final.drop(columns = ["Symbol","gc_content","gene_length", "Transcript count"],inplace=True)

    ## output csv with samples as columns, rows of genes 
    final.T.to_csv("BeWISE_scores.csv")
    """

}

