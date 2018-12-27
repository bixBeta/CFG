#!/usr/bin/python

import pandas as pd
import os;

os.getcwd( )
os.chdir('/Users/bix/Desktop/subset.0021/')

pd.read_csv("cpGs.genes.csv")
barCode = pd.read_csv("barCode.to.ID.csv")

d = barCode('val1').to_dict()
betas = pd.read_csv("0021.geneSubset.csv")

# betas.index = betas['cpGs']

betas.set_index('cpGs')


#df.rename(columns=dict(zip(barCode["val1"], [""])))
betas.rename(columns=dict(zip(barCode['val1'], barCode['val2'])))

#df1.columns=[df2['val2']]
#betas.columns=[barCode['val2']]
