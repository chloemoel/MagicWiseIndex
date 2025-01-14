#!/usr/bin/env nextflow

process COMBINE {
    container "oras://community.wave.seqera.io/library/pip_numpy_pandas:b978164578f33ac0"
    
    publishDir params.outdir, mode: 'symlink'

    input:
        path genetic_score
        path bewise_score
        
    output:
        path "CombinedBurdenEstimate.csv", emit: csv

    script:
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np

        genetic_score = pd.read_csv(${genetic_score}, index_col=0)
        bewise_score = pd.read_csv("${bewise_score}, index_col=0)

        common_samples = bewise_score.merge(genetic_score, left_index=True, right_index=True)
        common_genes = bewise_score.T.merge(genetic_score.T, left_index=True, right_index=True)

        common_genetic = genetic_score.loc[common_samples.index, common_genes.index]
        common_bewise = bewise_score.loc[common_samples.index, common_genes.index]
        
        combined_score = common_genetic + common_bewise

        combined_score.to_csv("CombinedBurdenEstimate.csv")
        """
}