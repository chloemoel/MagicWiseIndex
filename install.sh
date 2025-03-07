#!/bin/bash
#create singularity/apptainer containers with datasets and tools needed
cd bin

#vep
singularity pull --name vep.sif docker://ensemblorg/ensembl-vep
singularity exec vep.sif INSTALL.pl -c ../data -a cfp -s homo_sapiens -y GRCh37 --PLUGINS CADD

#python tools
singularity build tools.sif tools.def