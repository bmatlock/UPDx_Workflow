#!bin/python

import argparse
import pandas as pd
import statistics



parser = argparse.ArgumentParser()
parser.add_argument('--proportions_txt', type=str, required=True, help='The proportions text file that needs to be analyzed')
parser.add_argument('--cutoff_multiplier', type=str, required=True, help='The cutoff multiplier text file that needs to be created')
args = parser.parse_args()

proportions_library = open(args.proportions_txt)
proportions_content = proportions_library.read()
proportion_list = proportions_content.split(sep=",")

proportion_df = pd.DataFrame(proportion_list, columns=['Proportion of contam Reads'])

proportion_df['Proportion of contam Reads'] = pd.to_numeric(proportion_df['Proportion of contam Reads'])


average = proportion_df['Proportion of contam Reads'].mean()

SD = proportion_df['Proportion of contam Reads'].std(ddof=0)

print(SD)
print(average)

cutoff_multiplier = str(average + (4*SD))

with open (args.cutoff_multiplier, "w") as f:
    f.write(cutoff_multiplier)
