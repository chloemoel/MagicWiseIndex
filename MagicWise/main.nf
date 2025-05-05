#!/usr/bin/env nextflow

process CALCULATE_MAGICWISE {
    publishDir params.outdir, mode: 'copy'

    input:
        path genetic_score
        path bewise_score
        
    output:
        path "MagicWise.csv", emit: csv

    script:
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np

        genetic_score = pd.read_csv("${genetic_score}", index_col=0)
        bewise_score = pd.read_csv("${bewise_score}", index_col=0)

        common_samples = bewise_score.merge(genetic_score, left_index=True, right_index=True)
        common_genes = bewise_score.T.merge(genetic_score.T, left_index=True, right_index=True)

        common_genetic = genetic_score.loc[common_samples.index, common_genes.index]
        common_bewise = bewise_score.loc[common_samples.index, common_genes.index]

        phred_bewise = common_bewise.rank(method="dense", ascending=False).apply(lambda x: -10*(np.log10(1/x)))
        phred_genetic = common_genetic.rank(method="dense", ascending=False).apply(lambda x: -10*(np.log10(1/x)))
        
        combined_score = phred_genetic + phred_bewise

        combined_score.to_csv("MagicWise.csv")
        """
}
