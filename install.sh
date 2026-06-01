#!/bin/bash
cd bin/data

wget https://github.com/chloemoel/MagicWiseIndex/releases/download/v.0.0.2/dbNSFP4.9a.MagicWise.txt.gz
echo "Unzipping annotation info..."
tar -xvf dbNSFP4.9a.MagicWise.txt.gz

wget https://github.com/chloemoel/MagicWiseIndex/releases/download/v.0.0.2/probe_info.csv.gz
echo "Unzipping methylation probe info..."
tar -xvf probe_info.csv.gz

echo "All done!"
