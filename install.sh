#!/bin/bash
cd bin/data

wget https://github.com/chloemoel/MagicWiseIndex/releases/download/v.0.0.2/dbNSFP4.9a.MagicWise.txt.gz
echo "Unzipping annotation info..."
gunzip dbNSFP4.9a.MagicWise.txt.gz

wget https://github.com/chloemoel/MagicWiseIndex/releases/download/v.0.0.2/probe_info.csv.gz
echo "Unzipping methylation probe info..."
gunzip probe_info.csv.gz

echo "All done!"
