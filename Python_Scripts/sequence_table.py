#!/bin/python

from ast import arg
import string
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--Complete_Summary_Table', type=argparse.FileType('r'), required=True, help='The complete summary table that needs to be analyzed')
parser.add_argument('--Sequence_Summary_Table', type=argparse.FileType('a'), required=True, help='The new summary table that needs to be analyzed')
parser.add_argument('--Specimen_Name', type=str, required=True, help="The name of the specimen")
#parser.add_argument('--fasta_table', type=argparse.FileType('a'), required=True, help='The fasta table that needs to be created')
args = parser.parse_args()

Complete_df = pd.read_csv(args.Complete_Summary_Table)

temp1_fasta_table = Complete_df.drop(['% Assembled', 'E-Value', 'Cutoff Threshold', 'BLAST_Sequence','Cluster #', 'Total # of Sequences in Cluster', 'Total # of Sequences in Sample', 'PIDENT'], axis=1)
temp2_fasta_table = temp1_fasta_table.loc[temp1_fasta_table['Positive or Negative'] == 'Positive']
temp3_fasta_table= temp2_fasta_table.drop(['Positive or Negative', 'SseqID', 'Sequence'], axis=1)
sequences = temp2_fasta_table['SeqID']
fasta_table = pd.DataFrame({args.Specimen_Name:sequences}).reset_index(drop=True)
#SeqID = temp2_fasta_table['SeqID']
Sequence = temp2_fasta_table['Sequence']
SseqID = temp2_fasta_table['SseqID']
#fasta_table['SeqID'] = SeqID
fasta_table['Sequence'] = Sequence
fasta_table['SseqID'] = SseqID


fasta_table.to_csv(args.Sequence_Summary_Table, index=False)