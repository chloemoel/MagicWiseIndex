#!/usr/bin/env nextflow

process ANNOTATE_VCF {
    // pulls container from dockerhub
    container "ensemblorg/ensembl-vep"

    publishDir params.outdir, mode: 'symlink'

    input:
       path gvcf
       path genome_fasta
       path CADD
    
    output:
        path "annotated.merged.vcf", emit: vcf

    script:
    """
    vep \
        -a GRCh37 \
        --cache \
        --input_file ${gvcf} \
        --output_file "annotated.merged.vcf" \
        --fasta ${genome_fasta} \
        --dir_plugins ${CADD} \
        --plugin CADD \
        --vcf \
        --force_overwrite
    """
}

