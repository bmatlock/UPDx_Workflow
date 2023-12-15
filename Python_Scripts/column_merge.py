#!/bin/python

import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--count_summary', type=argparse.FileType('r'), help='The current count dataframe that needs to be merged')
parser.add_argument('--rep_seq_summary', type=argparse.FileType('r'), help='The current seqid dataframe that needs to be merged')
parser.add_argument('--full_summary', type=argparse.FileType('w'), help='The path to and name of overall summary table being created')
args = parser.parse_args()



count = pd.read_csv(args.count_summary)
seqid = pd.read_csv(args.rep_seq_summary)




test = pd.merge(seqid, count)


#ADDED THIS PART IN ORDER TO RESET THE INDEX AND ADD NEW HEADERS TO THE SUMMARY TABLE

#transpose, reset_index, transpose back, reset_index again
summary = test.T.reset_index().T.reset_index(drop=True)
#Create new column headers for your table
summary.columns = ['Cluster #', 'Nucleotides', 'SeqID', 'Total # of Sequences in Cluster']
#Drop the nucleotides column
summary.drop(columns=['Nucleotides'], inplace=True)

#Sum the total number of sequences for the sample and add it in a new column
summary['Total # of Sequences in Cluster'] = summary['Total # of Sequences in Cluster'].astype(int)

cluster_sum = summary['Total # of Sequences in Cluster'].sum()

summary['Total # of Sequences in Sample'] = cluster_sum


for ind,row in summary.iterrows():
    summary.loc[ind, "% Assembled"] = (row['Total # of Sequences in Cluster']/cluster_sum) * 100


summary.to_csv(args.full_summary, index=False)