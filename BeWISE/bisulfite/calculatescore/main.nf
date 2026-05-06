#!/usr/bin/env nextflow

process CALCULATE_BEWISE_BIS {
    publishDir params.outdir, mode: 'copy'
    
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

    probe_info = pd.read_csv("${probe_info}")
    m_values = pd.read_csv("${processed_m_values}")
    
    sample_names = m_values.columns[3:]
    score = m_values.merge(probe_info, left_on=["chr","start","stop"], right_on=["chr","start","stop"])
    print(score)
    score[sample_names] = score[sample_names].multiply(score["weights"],axis=0)


    # This collapses the cpgs (that are weighted by chromatin state) into a single value per gene
    gene_variables = ["Gene","gene_length","Transcript count"]
    collapse = score[list(sample_names) + gene_variables].groupby(["Gene","gene_length","Transcript count"]).agg('sum')

    samples = collapse[sample_names]
    length_row = collapse.index.get_level_values("gene_length")
    transcript_count_row = collapse.index.get_level_values("Transcript count")

    length = np.sqrt(length_row)
    transcript_count = np.sqrt(transcript_count_row)
    
    weighted_by_gene = samples.div(length, axis = "index").div(transcript_count, axis = "index")
    weighted_by_gene.index = collapse.index.get_level_values("Gene")

    ## output csv with samples as rows, genes as columns
    weighted_by_gene.T.to_csv("BeWISE_scores.csv")
    """

}
