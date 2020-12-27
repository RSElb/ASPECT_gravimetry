import numpy as np
from scipy.special import lpmv
import pandas as pd
import pyshtools as pysh
import matplotlib.pyplot as plt

inputfile = 'geoid_grid0.txt' # File to read the grid data from
outputfile = 'coeffs.txt' # File to write the SPH coefficient to

print('Reading in grid data from "%s"...' % inputfile)
df = pd.read_csv(inputfile, delim_whitespace=True, names=['Value'], dtype=np.float64)

values = np.array(df['Value'])
values = values.reshape((181, 361))

grid = pysh.SHGrid.from_array(values)

#fix, ax = grid.plot()
#plt.savefig('test.png', dpi=300)

print('Saving SPH coefficients in "%s"...' % outputfile)
coeffs = grid.expand(normalization='Ortho')
coeffs.to_file(outputfile)
