#!/usr/bin/env nextflow

process CALCULATE_BEMAGIC {
    publishDir params.outdir, mode: 'copy'
    
    input:
        path vcf
        path probe_info
        
    output:
        path "BeMAGIC_score.csv", emit: csv

    script:
    """
    #!/usr/bin/env python

    # Dependencies
    import pandas as pd
    import numpy as np
    from scipy.stats import *

    pd.set_option('future.no_silent_downcasting', True)
    
    # Load gene information
    gene_info = pd.read_csv("${probe_info}", usecols=["Gene", "gene_length", "Transcript count"])
    gene_info.drop_duplicates(inplace=True)    
    
    chunksize=10**6
    
    # Load the VCF file
    for vcf in pd.read_csv('${vcf}', sep='\\t', low_memory=False, chunksize=chunksize):

        # Filter based on the chromosome and ensure the index aligns
        vcf = vcf[~vcf["CHROM"].isin(["MT", "Y", "X"])]
  
        # Add CADD_RAW and Gene as its own column and explode any entries with multiple genes
        vcf["CADD_RAW"] = vcf["INFO"].apply(lambda x: x.split('|')[-1].split("CADD_raw_hg19=")[-1])
        vcf["CADD_RAW"] = vcf["CADD_RAW"].replace(to_replace='.', value = 0)
        vcf["Gene"] = vcf["INFO"].apply(lambda x: x.split('|')[0].split("genename=")[-1])
        vcf["Gene"] = vcf['Gene'].str.split(';').apply(lambda genes: list(set(genes)))
        vcf = vcf.explode("Gene", ignore_index=True)

        # Get the sample names (from the 9th column onward)
        sample_names = vcf.columns[~vcf.columns.isin(["CADD_RAW","Gene","INFO","CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT"])]  # These should be the sample names in the VCF file

        # Process the genotypes for each sample
        for col in sample_names:
            vcf[col] = vcf[col].apply(lambda x: x.split(":")[0])

        # Replace genotypes with numeric values
        vcf.replace(to_replace=['0/0'], value=0, inplace=True)
        vcf.replace(to_replace=['0/1'], value=1, inplace=True)
        vcf.replace(to_replace=['1/1'], value=2, inplace=True)
        vcf.replace(to_replace=['0/2'], value=None, inplace=True)
        vcf.replace(to_replace=['./.'], value=None, inplace=True)  

        # Drop rows with too many missing values
        thresh = int(len(sample_names) * 0.1)
        vcf.dropna(thresh=thresh, inplace=True)

        # Fill missing values with 0
        vcf = vcf.replace(to_replace=['', ' ', np.nan], value = 0)

        # Apply CADD score
        scored = vcf[sample_names].astype(int).mul(vcf["CADD_RAW"].astype(float),axis=0)
        scored = scored.set_index(vcf["Gene"])

        # Aggregate by gene
        grouped_gene = scored.groupby(["Gene"]).agg('sum')


        # Merge with gene info
        grouped_gene_info = grouped_gene.merge(gene_info, left_index=True, right_on="Gene")

        # Standardize gene length and transcript count
        length = zscore(grouped_gene_info["gene_length"])
        transcript = zscore(grouped_gene_info["Transcript count"])

        # Normalize the score by gene length and transcript count
        scored_gene = grouped_gene_info.drop(columns=gene_info.columns).T.astype(float) / transcript / length
        scored_gene.columns = grouped_gene_info["Gene"]

        # Save the final result to a CSV with genes as columns and samples as rows
        scored_gene.to_csv("BeMAGIC_score.csv", mode="a")
    """
}
