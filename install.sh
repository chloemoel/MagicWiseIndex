#!/bin/bash
cd bin/data

wget https://m-42bb44.6ba50.0ec8.data.globus.org/Shared/darbro-cytogenetics/Chloe/CombinedBurdenEstimate/bin/data/dbNSFP4.9a.MagicWise.txt.gz

echo "Unzipping annotation info..."
gunzip dbNSFP4.9a.MagicWise.txt.gz

echo "All done!"
