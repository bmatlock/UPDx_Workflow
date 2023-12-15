#!/bin/python

from ast import arg
import string
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cutoff_multiplier', type=str, required=True, help='The cutoff multiplier that needs to be used')
parser.add_argument('--Complete_Summary_Table', type=argparse.FileType('r'), required=True, help='The complete summary table that needs to be analyzed')
parser.add_argument('--Cutoff_Summary_Table', type=argparse.FileType('a'), required=True, help='The new summary table that needs to be analyzed')
parser.add_argument('--Sample_Parasite_Table', type=argparse.FileType('a'), required=True, help='The final parasite table that needs to be created for each table')
parser.add_argument('--markdown_summary_table', type=argparse.FileType('a'), required=True, help="The new r markdown output table that needs to be created")
parser.add_argument('--threshold_summary_text', type=str, required=True, help="The temporary threshold text file that needs to be creaeated for r markdown report")
parser.add_argument('--results_table', type=argparse.FileType('a'), required=True, help="The results tablee that needs to be created for r mardown report")
parser.add_argument('--specimen_name', type=str, required=True, help='The name of the current specimen being analyzed')
#parser.add_argument('--fasta_table', type=argparse.FileType('a'), required=True, help='The fasta table that needs to be created')
args = parser.parse_args()

Complete_df = pd.read_csv(args.Complete_Summary_Table, index_col=0)

cutoff_library = open(args.cutoff_multiplier)
cutoff_contents = cutoff_library.read()
cutoff_int = float(cutoff_contents)

Threshold_column = Complete_df['Total # of Sequences in Sample'].apply(lambda x: x * cutoff_int)

Complete_df['Cutoff Threshold'] = Threshold_column

Complete_df.reset_index(inplace=True, drop=True)

Complete_df['Cutoff Threshold'].fillna(20, inplace=True)

threshold = Complete_df.at[1, 'Cutoff Threshold'] 

threshold_int = float(threshold)

Complete_df.loc[Complete_df['Total # of Sequences in Cluster'] > threshold_int, 'Positive or Negative'] = 'Positive'
Complete_df.loc[Complete_df['Total # of Sequences in Cluster'] < threshold_int, 'Positive or Negative'] = 'Negative'

summary_df = Complete_df[Complete_df['Positive or Negative'] == 'Positive']


Final_summary_df = pd.DataFrame().assign(Parasite=summary_df['SseqID'])

parasite_hits = Final_summary_df['Parasite']
 
parasite_results_df=pd.DataFrame({args.specimen_name:parasite_hits}).reset_index(drop=True)

parasite_results_df.to_csv(args.Sample_Parasite_Table, index=False)

Complete_df.to_csv(args.Cutoff_Summary_Table, index=False)

markdown_summary_table = Complete_df.drop(['% Assembled', 'E-Value', 'Cutoff Threshold', 'BLAST_Sequence'], axis=1)
markdown_summary_table['Sequence'] = markdown_summary_table['Sequence'].str.wrap(15)

results_table = markdown_summary_table.drop(['Cluster #', 'Total # of Sequences in Cluster', 'Total # of Sequences in Sample', 'Sequence', 'PIDENT'], axis=1)

temp1_fasta_table = Complete_df.drop(['% Assembled', 'E-Value', 'Cutoff Threshold', 'BLAST_Sequence','Cluster #', 'Total # of Sequences in Cluster', 'Total # of Sequences in Sample', 'PIDENT', 'SseqID'], axis=1)
temp1_fasta_table['Sequence'] = temp1_fasta_table['Sequence'].str.wrap(500)

#temp2_fasta_table = temp1_fasta_table.loc[temp1_fasta_table['Positive or Negative'] == 'Positive']

#fasta_table = temp2_fasta_table.drop(['Positive or Negative'], axis=1)

threshold_str = str(threshold_int)

with open(args.threshold_summary_text, 'w') as f:
    f.write(threshold_str)

markdown_summary_table.to_csv(args.markdown_summary_table, index=False)
results_table.to_csv(args.results_table, index=False)

#fasta_table.to_csv(args.fasta_table, index=False, header=False)


