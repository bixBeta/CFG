#!/usr/bin/python

import pandas as pd
import os
import sys

os.getcwd( )
# os.chdir('/Users/bix/Desktop/subset.0021/') # change dir if needed


barCode = pd.read_csv('key.csv') # read in the key value combinations

betas = pd.read_csv('EPIC.csv') # read in the final dropped beta values with centrix id's as colnames

betas.set_index('samples') # specify string from the -2 argument file to use as index[i]

# rename based on specified key value relationship (manually change)
finalSet = betas.rename(columns=dict(zip(barCode['centrix'], barCode['Sample_Name']))) 

