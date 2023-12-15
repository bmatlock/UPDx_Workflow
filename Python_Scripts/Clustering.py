#!/bin/python

#!/bin/python
#blast_output /Users/brycematlock/Documents/UPDx_Workflow/Results/Experiment_10-05-23-15_49_20/FINAL_run_summary_10052023-15_49_20/Samples/2017_1742_2_S4_L001..summary/Tree_Building/unclear_hits/Clustering/BLAST/Results/M05039:164:000000000-D6HCM:1:1101:19570:9949.results.modified
#unclear_sequence_BLAST_summary_table /Users/brycematlock/Documents/unclear_sequence_BLAST_summary.csv
#unclear_sequence_fasta_table /Users/brycematlock/Documents/unclear_sequence_fasta_table.csv

import re
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('--blast_output', type=str, required=True, help='The full reference blast output that needs to be analyzed')
parser.add_argument('--unclear_sequence_BLAST_summary_table', type=argparse.FileType('a'), required=True, help='The new parasite BLAST output table that needs to be created')
args = parser.parse_args()

SeqID_Library = open(args.blast_output)
SeqID_Contents = SeqID_Library.read()
pattern = re.compile(r'>([A-Z0-9:-]+)\t[0-9]+\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+.*[0-9a-z-.]+')
matches = pattern.findall(SeqID_Contents)
seqid_list = []
for match in matches:
    seqid_list.append(match)
    SeqID_DF = pd.DataFrame(seqid_list, columns=['SeqID'])
    
Bitscore_Library = open(args.blast_output)
Bitscore_Contents = Bitscore_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t([0-9]+)\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+.*[0-9a-z-.]+')
matches = pattern.findall(Bitscore_Contents)
Bitscore_list = []
for match in matches:
    Bitscore_list.append(match)
    Bitscore_DF = pd.DataFrame(Bitscore_list, columns=['Bitscore'])   

Temporary_Table0 = pd.concat([SeqID_DF.reset_index(drop=True), Bitscore_DF.reset_index(drop=True)], axis=1)

Pident_Library = open(args.blast_output)
Pident_Contents = Pident_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9]+\t([0-9.]+)\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+.*[0-9a-z-.]+')
matches = pattern.findall(Pident_Contents)
pident_list = []
for match in matches:
    pident_list.append(match)
    Pident_DF = pd.DataFrame(pident_list, columns=['PIDENT'])
    
    

Temporary_Table1 = pd.concat([Temporary_Table0.reset_index(drop=True), Pident_DF.reset_index(drop=True)], axis=1)

SseqID_Library = open(args.blast_output)
SseqID_Contents = SseqID_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9]+\t[0-9.]+\t([A-Za-z0-9-_.():{}[\]<>*,]+)\t[A-Z-]+.*[0-9a-z-.]+')
matches = pattern.findall(SseqID_Contents)
sseqid_list = []
for match in matches:
    sseqid_list.append(match)
    SseqID_DF = pd.DataFrame(sseqid_list, columns=['SseqID'])

Temporary_Table2 = pd.concat([Temporary_Table1.reset_index(drop=True), SseqID_DF.reset_index(drop=True)], axis=1)


Sequence_Library = open(args.blast_output)
Sequence_Contents = Sequence_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9]+\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t([A-Z-]+).*[0-9a-z-.]+')
matches = pattern.findall(Sequence_Contents)
sequence_list = []
for match in matches:
    sequence_list.append(match)
    Sequence_DF = pd.DataFrame(sequence_list, columns=['BLAST_Sequence'])

Temporary_Table3 = pd.concat([Temporary_Table2.reset_index(drop=True), Sequence_DF.reset_index(drop=True)], axis=1)

Evalue_Library = open(args.blast_output)
Evalue_Contents = Evalue_Library.read()
pattern = re.compile(r'>[A-Z0-9:-]+\t[0-9]+\t[0-9.]+\t[A-Za-z0-9-_.():{}[\]<>*,]+\t[A-Z-]+\t([0-9a-z-.]+)')
matches = pattern.findall(Evalue_Contents)
evalue_list = []
for match in matches:
    evalue_list.append(match)
    Evalue_DF = pd.DataFrame(evalue_list, columns=['E-Value'])

unclear_sequence_BLAST_Summary_Table_temp1 = pd.concat([Temporary_Table3.reset_index(drop=True), Evalue_DF.reset_index(drop=True)], axis=1)
unclear_sequence_BLAST_Summary_Table_temp2 = unclear_sequence_BLAST_Summary_Table_temp1.sort_values(['Bitscore'], ascending=False)
unclear_sequence_BLAST_Summary_Table = unclear_sequence_BLAST_Summary_Table_temp2.head(40)
unclear_sequence_fasta_table = unclear_sequence_BLAST_Summary_Table.drop(['SeqID', 'Bitscore', 'PIDENT', 'E-Value', 'BLAST_Sequence'], axis=1)
unclear_sequence_BLAST_Summary_Table.to_csv(args.unclear_sequence_BLAST_summary_table)