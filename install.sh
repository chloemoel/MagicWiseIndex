#!/bin/bash
cd bin/

# Check if the user passed "docker" as an argument
if [ "$1" == "docker" ]; then
    echo "Installing tools with Docker..."

    # Build Docker image for Python tools (based on tools.def)
    docker build -f tools_docker.def -t tools_docker .
    
else
    echo "Installing tools with Singularity..."
    # Build Python tools container
    singularity build tools.sif tools.def
fi
