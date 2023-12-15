#!/bin/python

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--count_df', type=argparse.FileType('r'), required=True, help='The count summary dataframe that needs to be analyzed')
parser.add_argument('--clusters_filtered', type=str, required=True, help='The list of clusters that have >20 sequences in a text file')
args = parser.parse_args()

count_df = pd.read_csv(args.count_df)

count_df.columns = ['Cluster', 'Count']

filtered_df = count_df[count_df['Count'] >= 20]

with open (args.clusters_filtered, 'w') as f:
    clusters = filtered_df['Cluster'].to_string(header=False, index=False)
    f.write(clusters)

