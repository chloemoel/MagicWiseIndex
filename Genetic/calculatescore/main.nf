#!/usr/bin/env nextflow

process CALCULATE_GENE_SCORE {
    container "oras://community.wave.seqera.io/library/pip_numpy_pandas:b978164578f33ac0"
    
    publishDir params.outdir, mode: 'symlink'

    input:
        path vcf
        path probe_info
        path sample_sheet
        
    output:
        path("genetic_score.csv"), emit: csv

    script:
        """
        #!/usr/bin/env python

        # Dependencies
        import pandas as pd
        import numpy as np

        def get_vcf_names(vcf_path):
            with open(vcf_path, "rt") as ifile:
                for line in ifile:
                    if line.startswith("#CHROM"):
                        vcf_names = line.strip('#\n').split('\t')
                        break
            ifile.close()
            return vcf_names


        names = get_vcf_names('${vcf}')
        vcf = pd.read_csv('${vcf}', comment='#', sep='\s+', header=None, names=names, low_memory=False)

        
        def split_vep_fields(vep_col):
            fields = vep_col.split("|")
            vep_dict = {
                "Consequence": fields[1],
                "IMPACT": fields[2],
                "SYMBOL": fields[3],
                "Gene": fields[4],
                # "Feature_type": fields[5],
                # "Feature": fields[6],
                # "BIOTYPE": fields[7],
                # "EXON": fields[8],
                # "INTRON": fields[9],
                # "HGVSc": fields[10],
                # "HGVSp": fields[11],
                # "cDNA_position": fields[12],
                # "CDS_position": fields[13],
                # "Protein_position": fields[14],
                # "Amino_acids": fields[15],
                # "Codons": fields[16],
                # "Existing_variation": fields[17],
                # "DISTANCE": fields[18],
                # "STRAND": fields[19],
                # "FLAGS": fields[20],
                # "VARIANT_CLASS": fields[21],
                # "SYMBOL_SOURCE": fields[22],
                # "HGNC_ID": fields[23],
                # "ENSP": fields[24],
                "SIFT": fields[25],
                "PolyPhen": fields[26],
                # "HGVS_OFFSET": fields[27],
                # "MAX_AF": fields[28],
                # "MAX_AF_POPS": fields[29],
                # "CLIN_SIG": fields[30],
                # "SOMATIC": fields[31],
                # "PHENO": fields[32],
                "CADD_PHRED": fields[33],
                "CADD_RAW": fields.pop()
            }

            return vep_dict

        vep_dict = split_vep_fields(vcf["INFO"][0])
        cols = list(vep_dict.keys())
        vcf[cols] = vcf["INFO"].apply(lambda x: pd.Series(split_vep_fields(x)))


        # Filter based on the chromosome and ensure the index aligns
        vcf = vcf[~vcf["CHROM"].isin(["MT", "Y", "X"])]

        # Identify sample names and info columns
        sample_names_full = [col for col in vcf.columns if col.startswith('2')]
        info = [col for col in vcf.columns if not col.startswith('2')]

        # Create a list to hold Series for each sample's genotype
        genotype_series = []

        # Extract genotypes for each sample
        for col in sample_names_full:
            # Split the genotype information and take the first part
            genotype_series.append(vcf[col].apply(lambda x: x.split(":")[0]))

        # Concatenate all genotype Series into a DataFrame
        genotypes = pd.concat(genotype_series, axis=1)
        genotypes.index = vcf["ID"]

        # Replace blanks with 0
        genotypes = genotypes.replace(to_replace=['', ' ', None], value=0)

        # Generate new column names
        sample_names = [x.split(".")[0] for x in sample_names_full]
        columns = sample_names

        # Set new column names
        genotypes.columns = columns

        genotypes = genotypes.merge(vcf[info], left_index = True, right_on = "ID")
        genotypes.fillna(0,inplace=True)
        genotypes.replace({"":0," ":0}, inplace=True)
        genotypes = genotypes.assign(genetic_score =  genotypes["CADD_RAW"].astype(float))

        sample_group = pd.read_csv("${sample_sheet}")

        genotypes.replace(to_replace=['0/0'], value=0, inplace=True)
        genotypes.replace(to_replace=['0/1'], value=1, inplace=True)
        genotypes.replace(to_replace=['1/1'], value=2, inplace=True)
        genotypes.replace(to_replace=['0/2'], value=None, inplace=True)
        genotypes.replace(to_replace=['./.'], value=None, inplace=True)  

        thresh = int(len(sample_names) * 0.1)
        genotypes.dropna(thresh=thresh , inplace=True)

        #Fill missing values with 0
        ####should this maybe be median?####
        genotypes.fillna(0, axis=0, inplace=True) 

        scored = genotypes[sample_names].T.astype(int) * genotypes["genetic_score"]
        scored_snps = pd.concat([genotypes[["CHROM", "POS", "ID", "Gene", "SYMBOL"]].T, scored]).T 

        grouped_gene = scored_snps.groupby(["Gene","CHROM"]).agg('sum').reset_index()

        gene_info = pd.read_csv("${probe_info}, use_cols=["Symbol","gene_length","Transcript count"])
        gene_info.drop_duplicates(inplace=True)

        grouped_gene_info = grouped_gene.merge(gene_info,left_on="Gene",right_on="Symbol",how="left")

        length = zscore(grouped_gene_info["length"])
        transcript = zscore(grouped_gene_info["Transcript count"])

        scored_gene = grouped_gene.astype(float)) / transcript / length

        scored_gene.to_csv("genetic_score.csv")
        """
}

