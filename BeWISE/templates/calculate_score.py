#!/usr/bin/env python3

# Dependencies
import argparse
import pandas as pd
import numpy as np


def main():
 
    parser = argparse.ArgumentParser(description="Calculate BeWISE",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("probe_info", 
                        help='''csv file with probe information: 
                        (ID,Gene,gc_content,Symbol,gene_length,Transcript count,weights)''')
    parser.add_argument("m_values",
                        help="CLEANED csv file with sentrix IDs on top, and probe names on the side")

    args = parser.parse_args()
    config = vars(args)

    probe_info = pd.read_csv(config["probe_info"])
    m_values = pd.read_csv(config["m_values"])

    ## Merge m values with probe information. Do "right" merge to filter down probes
    score = m_values.merge(probe_info["ID", "weights"], left_index=True, right_on="ID", how="right")

    sample_columns = score.drop(columns=["ID","weights"])
    cpg_weights = score["weights"]
    idx = score["ID"]

    # Multiply cpg values with cpg_weights
    weighted_cpgs = sample_columns * cpg_weights
    weighted_cpgs.index = idx

    #this merges our weighted cpg site with the info on the gene
    full = weighted_cpgs.merge(probe_info, left_index=True, right_on="ID", how="right")
    full.drop(columns="weights", inplace=True)
    full.set_index("ID", inplace=True)
    full.drop_duplicates(inplace=True)

    # This collapses the cpgs (that are weighted by chromatin state) into a single value per gene
    collapse = full.groupby(["Symbol", "gc_content","gene_length", "Transcript count"]).agg('sum')

    ## output csv with samples as columns, rows of genes 
    collapse.to_csv("BeWISE_scores.csv", index=False)

if __name__ == "__main__":
	main()





