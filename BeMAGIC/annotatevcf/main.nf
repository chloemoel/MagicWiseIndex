#!/usr/bin/env nextflow

process ANNOTATE_VCF {
    publishDir params.outdir, mode: 'copy'

    input:
        path dbnsfp
        path vcf
    
    output:
        path "${vcf}_ANNOTATED", emit: vcf

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np
    import gzip

    # Load dbNSFP (assuming it fits in memory)
    dbnsfp = pd.read_csv('${dbnsfp}', sep='\\t', dtype={'CHROM': 'Int64', 'POS': 'Int64'})
    dbnsfp.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    
    # Find number of rows to skip
    row_skip = 0
    with gzip.open("${vcf}", "rt") as ifile:
        for line in ifile:
            if not line.startswith("#CHROM"):
                row_skip += 1
            else:
                break


    # Load the VCF file in chunks to handle large datasets
    chunksize = 10**6  # Adjust this based on available memory

    round_time = 1

    # Read VCF in chunks to avoid loading everything into memory at once
    for chunk in pd.read_csv('${vcf}', skiprows=row_skip, sep='\\t', header=0, low_memory=False, chunksize=chunksize, compression="gzip"):
        # Remove duplicated columns and clean up the VCF data
        chunk = chunk.loc[:, ~chunk.columns.duplicated()]  # Takes first instance of duplicated column
        chunk.drop(columns=['ID', 'QUAL', 'FILTER', 'INFO'], inplace=True)
        chunk.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
        
        # Apply the necessary transformations if needed (ie chr1 to 1)
        if chunk['CHROM'].dtype == object:
            chunk['CHROM'] = chunk['CHROM'].apply(lambda x: x.split('chr')[-1])
        
        # Keep only autosomes        
        chunk = chunk[chunk['CHROM'].isin([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"])]

        # Ensure chromosomes and positions  are integers for merging
        chunk['CHROM'] = chunk['CHROM'].astype('Int64')
        chunk['POS'] = chunk['POS'].astype('Int64')

        # Merge VCF chunk with dbNSFP (done per chunk)
        annotated_chunk = chunk.merge(dbnsfp, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')
        
        # Write out annotated data incrementally to avoid memory overload. only add header on first round
        if round_time == 1:
            round_time += 1
            annotated_chunk.to_csv("${vcf}_ANNOTATED", sep='\\t', mode='a',header=True, index=False)
        else:
            annotated_chunk.to_csv("${vcf}_ANNOTATED", sep='\\t', mode='a',header=False, index=False)
    """
    
}
