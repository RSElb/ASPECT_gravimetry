import numpy as np
import pandas as pd
import sys as sys

try:
    inputfile1 = str(sys.argv[1])
    inputfile2 = str(sys.argv[2])
    outputfile = str(sys.argv[3])
except:
    print('Type two different input grid files and one output file!')

print('Reading data from "%s"...' % inputfile1)
df1 = pd.read_csv(inputfile1, delim_whitespace=True, names=['Value'], dtype=np.float64)

print('Reading data from "%s"...' % inputfile2)
df2 = pd.read_csv(inputfile2, delim_whitespace=True, names=['Value'], dtype=np.float64)

values1 = np.array(df1['Value'], dtype=np.float64)
values2 = np.array(df2['Value'], dtype=np.float64)

print('Substracting data in "%s" from "%s"...' % (inputfile2, inputfile1))
output = values1 - values2

print('Saving output grid in "%s"...' % outputfile)
np.savetxt(outputfile, output, delimiter='\t', fmt='%.16e')
