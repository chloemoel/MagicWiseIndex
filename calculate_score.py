#!/usr/bin/env python3

# Dependencies
import argparse
import pandas as pd
import numpy as np


def main():
 
    parser = argparse.ArgumentParser(description="Calculate BeWISE",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("probe_info", 
                        help='''csv file with probe informatino 
                        (probe name, probe_weight, transcript count, gc content, length, gene_name)''')
    parser.add_argument("m_values",
                        help="CLEANED csv file with sentrix IDs on top, and probe names on the side")

    args = parser.parse_args()
    config = vars(args)

    probe_info = pd.read_csv(config["probe_info"])
    m_values = pd.read_csv(config["m_values"])

    ## Merge m values with probe information 
    score = m_values.merge(probe_info)

    ## TODO: multiply m value by probe weight

    ## Group weighted probes by gene
    score = score.groupby(gene_name,method=sum)

    ## divide gene level score by gene info
    score = 

    ## output csv with samples as columns, rows of genes 
    pd.write_csv("BeWISE_scores.csv")

if __name__ == "__main__":
	main()





