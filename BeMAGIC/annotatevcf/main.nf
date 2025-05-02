#!/usr/bin/env nextflow

process ANNOTATE_VCF {
    publishDir params.outdir, mode: 'copy'

    input:
        path dbnsfp
        path vcf
    
    output:
        path "annotated.merged.vcf", emit: vcf

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np

    # Number of rows to skip before encountering the header in VCF file
    row_skip = 0

    # Find the header line to determine the row_skip value
    with open("${vcf}", "rt") as ifile:
        for line in ifile:
            if not line.startswith("#CHROM"):
                row_skip += 1
            else:
                break

    # Load dbNSFP (assuming it fits in memory)
    dbnsfp = pd.read_csv('${dbnsfp}', sep='\\t', header=None, dtype={'CHROM': 'Int64', 'POS': 'Int64'})
    dbnsfp.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

    # Load the VCF file in chunks to handle large datasets
    chunksize = 10**6  # Adjust this based on available memory

    # Read VCF in chunks to avoid loading everything into memory at once
    for chunk in pd.read_csv('${vcf}', skiprows=row_skip, sep='\\t', header=0, low_memory=False, chunksize=chunksize):
        # Remove duplicated columns and clean up the VCF data
        chunk = chunk.loc[:, ~chunk.columns.duplicated()]  # Takes first instance of duplicated column
        chunk.drop(columns=['ID', 'QUAL', 'FILTER', 'INFO'], inplace=True)
        chunk.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
        # Filter out rows with undesired chromosomes
        chunk = chunk[~chunk['CHROM'].isin(['MT', 'Y', 'X'])]                                                                        
        # Apply the necessary transformations
        if chunk['CHROM'].dtype == object:
            chunk['CHROM'] = chunk['CHROM'].apply(lambda x: x.split('chr')[-1]).astype('Int64')
        chunk['POS'] = chunk['POS'].astype('Int64')

        # Merge VCF chunk with dbNSFP (done per chunk)
        annotated_chunk = chunk.merge(dbnsfp, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')
        # Write out annotated data incrementally to avoid memory overload
        annotated_chunk.to_csv("annotated.merged.vcf", sep='\\t', mode='a', header=False, index=False)
    """
    
}

