#!/bin/bash

# Check if the user passed "docker" as an argument
if [ "$1" == "docker" ]; then
    echo "Installing tools with Docker..."

    # Build Docker image for VEP
    docker pull ensemblorg/ensembl-vep
    docker run --rm -v "$(pwd)/data:/data" ensemblorg/ensembl-vep INSTALL.pl -c /data -a cfp -s homo_sapiens -y GRCh37 --PLUGINS CADD

    # Build Docker image for Python tools (based on tools.def)
    docker build -f tools.def -t tools .
else
    echo "Installing tools with Singularity..."
    # Pull and execute the VEP container
    singularity pull --name vep.sif docker://ensemblorg/ensembl-vep
    singularity exec vep.sif INSTALL.pl -c ../data -a cfp -s homo_sapiens -y GRCh37 --PLUGINS CADD
    
    # Build Python tools container
    singularity build tools.sif tools.def
fi
