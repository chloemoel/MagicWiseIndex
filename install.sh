#!/bin/bash
cd bin/data

wget https://m-42bb44.6ba50.0ec8.data.globus.org/Shared/darbro-cytogenetics/Chloe/CombinedBurdenEstimate/bin/data/dbNSFP4.9a.MagicWise.txt.gz

echo "Unzipping annotation info..."
gunzip dbNSFP4.9a.MagicWise.txt.gz


wget https://m-42bb44.6ba50.0ec8.data.globus.org/Shared/darbro-cytogenetics/Chloe/CombinedBurdenEstimate/bin/data/probe_info.csv.gz
echo "Unzipping methylation probe info..."
gunzip probe_info.csv.gz

echo "All done!"
