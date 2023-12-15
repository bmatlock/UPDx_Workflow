#!/bin/python
#Complete_Summary_Table /Users/brycematlock/Documents/UPDx_Workflow/FINAL_run_summary_11182022/Negatives/Ground_chicken_S25_L..summary/COMPLETE_Summary_Table.csv
#proportion_txt /Users/brycematlock/Documents/UPDx_Workflow/FINAL_run_summary_11182022/Negatives/Ground_chicken_S25_L..summary/proportion.txt


import argparse
import pandas as pd



parser = argparse.ArgumentParser()
parser.add_argument('--Complete_Summary_Table', type=argparse.FileType('r'), required=True, help='The full complete summary table that needs to be analyzed')
parser.add_argument('--proportion_txt', type=str, required=True, help='The text file that will hold the calculated proportion' )
parser.add_argument('--list_of_parasites',type=str, required=True, help='The text file that will contain the list of parasites used for cut-off calculation')
parser.add_argument('--parasite_table', type=argparse.FileType('a'), required=True, help='The table that will contain the list of parasites used for r-markdown report')
args = parser.parse_args()

pd.set_option('display.max_colwidth', None)

table = pd.read_csv(args.Complete_Summary_Table)

table.fillna('No parasites found', inplace=True)

Total = table.at[0, 'Total # of Sequences in Sample']

list = table[['SseqID', 'Sequence']]

list['Sequence'] = list['Sequence'].str.wrap(15)

list.fillna('No parasites found', inplace=True)

contam_reads = table['Total # of Sequences in Cluster'].sum()

final_contam_sum = (contam_reads)

proportion = str(final_contam_sum / Total)

with open (args.proportion_txt, "w") as f:
    f.write(proportion)

with open (args.list_of_parasites, "w") as f:
    f.write(str(list))

list.to_csv(args.parasite_table)

