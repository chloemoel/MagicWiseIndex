#!/usr/bin/env nextflow

process ANNOTATE_VCF {
    publishDir params.outdir, mode: 'symlink'
    container "bin/vep.sif"

    input:
        path cache_dir
        path vcf
        path genome_fasta
        path cadd
        path cadd_tbi
    
    output:
        path "annotated.merged.vcf", emit: vcf

    script:
    """
    vep \
    -a GRCh37 \
    --cache \
    --dir ${cache_dir} \
    --input_file ${vcf} \
    --output_file annotated.merged.vcf \
    --fasta ${genome_fasta}\
    --plugin CADD,${cadd} \
    --vcf \
    --force_overwrite
    """
}

