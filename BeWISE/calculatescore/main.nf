#!/usr/bin/env nextflow

process CALCULATE_BEWISE {
    publishDir params.outdir, mode: 'symlink'
    container "bin/tools.sif"
    
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
    from scipy.stats import *

    probe_info = pd.read_csv("${probe_info}", index_col=0)
    m_values = pd.read_csv("${processed_m_values}", index_col=0)

    ## Merge m values with probe information so weights are in correct order
    score = m_values.merge(probe_info[["weights"]], left_index=True, right_index=True)

    sample_columns = score.drop(columns=["weights"])
    cpg_weights = score["weights"]

    # Multiply cpg values with cpg_weights
    weighted_cpgs = sample_columns.multiply(cpg_weights, axis="index")

    #this merges our weighted cpg site with the info on the gene
    full = weighted_cpgs.merge(probe_info.loc[:, probe_info.columns != "weights"], left_index=True, right_index=True)
    full.drop_duplicates(inplace=True)

    # This collapses the cpgs (that are weighted by chromatin state) into a single value per gene
    collapse = full.groupby(["Gene","gene_length","Transcript count"]).agg('sum')

    samples = collapse.reset_index(drop=True)
    length = zscore(collapse.index.get_level_values("gene_length"))
    transcript_count = zscore(collapse.index.get_level_values("Transcript count"))

    weighted_by_gene = samples.div(length, axis = "index").div(transcript_count, axis = "index")
    weighted_by_gene.index = collapse.index.get_level_values("Gene")

    ## output csv with samples as columns, rows of genes 
    weighted_by_gene.T.to_csv("BeWISE_scores.csv")
    """

}

