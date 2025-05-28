#!/bin/bash
cd bin/

# Check if the user passed "docker" as an argument
if [ "$1" == "docker" ]; then
    echo "Installing tools with Docker..."

    # Build Docker image for Python tools (based on tools.def)
    docker build -f tools_docker.def -t tools_docker containers/
    
else
    echo "Installing tools with Singularity..."
    # Build Python tools container
    singularity build containers/tools.sif tools.def
fi

echo "Installing annotation info..."
cd data/
wget https://m-42bb44.6ba50.0ec8.data.globus.org/Shared/darbro-cytogenetics/Chloe/CombinedBurdenEstimate/bin/data/dbNSFP4.9a.MagicWise.txt.gz

echo "Unzipping annotation info..."
gunzip dbNSFP4.9a.MagicWise.txt.gz

echo "All done!"