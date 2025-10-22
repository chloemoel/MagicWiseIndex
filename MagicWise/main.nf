#!/usr/bin/env nextflow

process CALCULATE_MAGICWISE {
    publishDir params.outdir, mode: 'copy'

    input:
        path genetic_score
        path bewise_score
        path ens_to_gene 
    output:
        path "MagicWise_ens_ID.csv"
        path "MagicWise_genename.csv"

    script:
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np
        import csv

        bemagic = pd.read_csv("${genetic_score}", index_col=0)
        bemagic.index = bemagic.index.astype(str)
        bewise = pd.read_csv("${bewise_score}", index_col=0)
        bewise.index = bewise.index.astype(str)

        phred_bewise = bewise.rank(method="dense", ascending=False, axis=1).apply(lambda x: -10*(np.log10(x/len(bewise.columns))))
        phred_bemagic = bemagic.rank(method="dense", ascending=False, axis=1).apply(lambda x: -10*(np.log10(x/len(bemagic.columns))))

        common_samples = phred_bewise[phred_bewise.index.isin(phred_bemagic.index)].index
        common_genes = phred_bewise.T[phred_bewise.T.index.isin(phred_bemagic.T.index)].index

        combined_score = phred_bemagic.loc[common_samples,common_genes] + phred_bewise.loc[common_samples,common_genes]
        
        combined_score.to_csv("MagicWise_ens_ID.csv")
        
        data_dict = {}

        with open("${ens_to_gene}", mode="r", newline="") as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header row (remove if no header)

            for row in reader:
                key = row[0]
                value = row[1]
                data_dict[key] = value
        combined_score.columns = combined_score.columns.map(data_dict)
        combined_score.to_csv("MagicWise_genename.csv")	
		
        """
}
