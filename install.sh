#!/bin/bash
home_dir=$(pwd)
cd bin/

# Check if the user passed "docker" as an argument
if [ "$1" == "docker" ]; then
    echo "Installing tools with Docker..."

    # Build Docker image for VEP
    docker pull ensemblorg/ensembl-vep

    docker run -t -i -v $home_dir/data/vep_data:/data ensemblorg/ensembl-vep INSTALL.pl -c $home_dir/data/vep_data -a cfp -s homo_sapiens -y GRCh37 --PLUGINS CADD

    # Build Docker image for Python tools (based on tools.def)
    docker build -f tools_docker.def -t tools .
else
    echo "Installing tools with Singularity..."
    # Pull and execute the VEP container
    singularity pull --name vep.sif docker://ensemblorg/ensembl-vep
    singularity exec vep.sif INSTALL.pl -c $home_dir/data/vep_data -a cfp -s homo_sapiens -y GRCh37 --PLUGINS CADD
    
    # Build Python tools container
    singularity build tools.sif tools.def
fi