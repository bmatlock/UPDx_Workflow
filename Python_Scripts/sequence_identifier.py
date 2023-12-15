#!/bin/python


import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--fasta_file', type=str, required=True, help='The fasta file that needs to be analyzed')
parser.add_argument('--summary_table', type=argparse.FileType('r'), required=True, help='The Full Summary Table that needs to be updated')
parser.add_argument('--new_summary_table', type=argparse.FileType('a'), required=True, help='The New summary Table that needs to be created')
args = parser.parse_args()


SeqID_Library = open(args.fasta_file)
SeqID_Library_Contents = SeqID_Library.read()
seqID_pattern = re.compile(r'>([A-Z0-9:-]+)\s[0-9A-Z:]+\n')
matches = seqID_pattern.findall(SeqID_Library_Contents)
SeqID_list = []
for match in matches:
    SeqID_list.append(match)
    SeqID_DF = pd.DataFrame(SeqID_list, columns=['SeqID'])




Sequence_Library = open(args.fasta_file)
Sequence_Library_Contents = Sequence_Library.read()
sequence_pattern = re.compile(r'\n([A-Z]+)')
matches = sequence_pattern.findall(Sequence_Library_Contents)
sequence_list = []
for match in matches:
    sequence_list.append(match)
    sequence_DF = pd.DataFrame(sequence_list, columns=['Sequence'])
    
Sequence_summary_table = pd.concat([SeqID_DF.reset_index(drop=True), sequence_DF.reset_index(drop=True)], axis=1,)

#Sequence_summary_table.to_csv(args.sequence_table, index=False)

summary_DF = pd.read_csv(args.summary_table, sep=',')
full_summary_DF = pd.merge(summary_DF, Sequence_summary_table, on="SeqID")
full_summary_DF['Sequence'] = full_summary_DF['Sequence'].str.wrap(15)
full_summary_DF.sort_values(by="Total # of Sequences in Cluster", ascending=False, axis='index', inplace=True)
full_summary_DF.to_csv(args.new_summary_table, index=False)






