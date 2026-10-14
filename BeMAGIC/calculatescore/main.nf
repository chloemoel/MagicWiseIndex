#!/usr/bin/env nextflow

process CALCULATE_BEMAGIC {
    publishDir params.outdir, mode: 'copy'

    input:
    path vcf_files_nf_obj
    path probe_info

    output:
    path "BeMAGIC_score.csv"

    script:
    def vcf_files = vcf_files_nf_obj.collect { f -> "\"${f}\"" }.join(", ")
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

    gene_dict = {}

    # Read through all vcfs
    for vcf_path in [${vcf_files}]:

        # Check for duplicate sample names using the raw header line
        with open(vcf_path) as fh:
            raw_cols = fh.readline().rstrip("\\n").split("\\t")
        if len(raw_cols) != len(set(raw_cols)):
            print("WARNING: You have duplicate sample names! Only using first instance")

        header = pd.read_csv(vcf_path, sep='\\t', nrows=1, header=0)

        chunksize = 10**6

        # Load the VCF file in chunks
        for vcf in pd.read_csv(vcf_path, sep='\\t', low_memory=False, chunksize=chunksize, usecols=header.columns):

            # Filter based on the chromosome
            vcf = vcf[~vcf["CHROM"].isin(["MT", "Y", "X"])]

            # Add CADD_RAW and Gene as their own columns and explode entries with multiple genes
            vcf["CADD_RAW"] = vcf["INFO"].apply(lambda x: x.split('|')[-1].split("CADD_raw_hg19=")[-1])
            vcf["CADD_RAW"] = vcf["CADD_RAW"].replace(to_replace='.', value=0)
            vcf["Gene"] = vcf["INFO"].apply(lambda x: x.split('|')[1].split("ensembl_geneid=")[-1])
            vcf["Gene"] = vcf['Gene'].str.split(';').apply(lambda genes: list(set(genes)))
            vcf = vcf.explode("Gene", ignore_index=True)

            # Sample names
            sample_names = vcf.columns[~vcf.columns.isin(["CADD_RAW","Gene","INFO","CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT"])]

            # Keep only the genotype portion of each sample
            for col in sample_names:
                vcf[col] = vcf[col].apply(lambda x: x.split(":")[0])

            # Replace genotypes with numeric values
            vcf.replace(to_replace=r"0[/|]0", value=0, regex=True, inplace=True)
            vcf.replace(to_replace=r"0[/|]1", value=1, regex=True, inplace=True)
            vcf.replace(to_replace=r"1[/|]0", value=1, regex=True, inplace=True)
            vcf.replace(to_replace=r"1[/|]1", value=2, regex=True, inplace=True)
            vcf.replace(to_replace=r"0[/|]2", value=None, regex=True, inplace=True)
            vcf.replace(to_replace=r"[.][/|][.]", value=None, regex=True, inplace=True)

            # Drop rows with too many missing values
            thresh = int(len(sample_names) * 0.1)
            vcf.dropna(thresh=thresh, inplace=True)

            # Fill missing values with 0
            vcf = vcf.replace(to_replace=['', ' ', np.nan], value=0)

            # Apply CADD score
            scored = vcf[sample_names].astype(int).mul(vcf["CADD_RAW"].astype(float), axis=0)
            scored = scored.set_index(vcf["Gene"])

            # Aggregate by gene for this chunk
            grouped_gene = scored.groupby(["Gene"]).agg('sum')
            grouped_gene_info = grouped_gene.merge(gene_info, left_index=True, right_on="Gene")

            # Standardize gene length and transcript count
            length = np.sqrt(grouped_gene_info["gene_length"])
            transcript = np.sqrt(grouped_gene_info["Transcript count"])

            # Normalize the score by gene length and transcript count
            scored_gene = grouped_gene_info.drop(columns=gene_info.columns).T.astype(float) / transcript / length
            scored_gene.columns = grouped_gene_info["Gene"]

            gene_dict.update(scored_gene.to_dict("series"))

    # Write once: genes as columns, samples as rows
    print(gene_dict) 
    final_gene_df = pd.DataFrame.from_dict(gene_dict).fillna(0)
    final_gene_df.to_csv("BeMAGIC_score.csv")
    """
}
