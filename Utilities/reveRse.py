# python utility to create  reverse complements for NextSeq 500's i5 adapters
# written by Faraz Ahmed 
# FA180921

# version 1.1 

import os
import sys
from Bio.Seq import Seq


filename = sys.argv[-1]
sys.stdout = open('rcIndices', 'w')


file = open(filename, "r")                  # open file handle, read only
index = []                                  # initializes a list
for line in file:
    line = line.rstrip()                    # remove /n
    index.append( line )                    # add each line as a list element

file.close()                                # close file handle

for x in range(len(index)):
    dna = Seq(index[x])
    print (dna.reverse_complement())


os.rename('rcIndices', 'ReverseComplement_' + filename)






