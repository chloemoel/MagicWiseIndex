#!/usr/bin/env python3

# Dependencies
import argparse
import pandas as pd
import numpy as np
from inmoose.pycombat import pycombat_norm

def main():
 
    parser = argparse.ArgumentParser(description="Process methylation data for BeWISE",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-a,--additional_file", help="additional sample data for batch effects")
    parser.add_argument("-b","--batch", help = "List of columns from either file to use for batch correction")


    parser.add_argument("sample_sheet", 
                        help='''Assumes user has an Illumina style sample sheet in the format:
                        Sample_Name,Sample_Well,Sample_Plate,Sample_Group,Pool_ID,Sentrix_ID,Sentrix_Position''')
    parser.add_argument("m_values",
                        help=" a csv file with sentrix IDs on top, and probe names on the side")

    args = parser.parse_args()
    config = vars(args)

    m_values = pd.read_csv(config["m_values"])

    ## Check to see if the first lines are from the sample sheet or if already cleaned.
    with open(config["sample_sheet"], "r") as s:
        if s.readline().split(",") == " ":
            sample_info = pd.read_csv(config["sample_sheet"], skiprows=7)
        else:
            sample_info = pd.read_csv(config["sample_sheet"])

    # filter out sites with more than 10% of samples missing that probe. then fill NA with probe mean
    thresh = int(len(m_values) * 0.1)
    m_values.dropna(thresh=thresh , inplace=True, axis = 1)

    #fills NA values with the mean of the probe
    m_values.fillna(m_values.mean(), axis=0, inplace=True)

    if config["batch"].isnull() == False:
        if config["additional_file"].isnull == False:
            additional_data = pd.read_csv(config["additional_file"])
            m_values = m_values.merge(additional_data)

        batch_columns = np.array(config["batch"])

        #perform batch correction
        for b in batch_columns:
            ##TODO: fix m value dataframe so its just m values
            m_values = pycombat_norm(m_values.drop(columns=batch_columns), m_values[b])
            print("Batch correction for", b ,"done")


    pd.write_csv(m_values,"m_values_processed.csv")

if __name__ == "__main__":
	main()