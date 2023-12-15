#!/bin/python

import pandas as pd
import os
from glob import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--parasite_hit_tables', type=str, required=True, help='The parasite hit tables for each sample that needs to be combined')
parser.add_argument('--parasite_summary_file', type=argparse.FileType('a'), required=True, help='The final parasite summary output thtat needs to be created')
args = parser.parse_args()


summary_tables = [os.path.basename(x) for x in glob(os.path.join(args.parasite_hit_tables, '*.csv'))]
merged_summary_tables = pd.concat(pd.read_csv(summary_table).assign(sample = summary_table)
for summary_table in summary_tables)

merged_summary_tables.to_csv(args.parasite_summary_file)