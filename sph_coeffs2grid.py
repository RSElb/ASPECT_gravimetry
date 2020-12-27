import numpy as np
from scipy.special import lpmv
import pandas as pd
import pyshtools as pysh

# degree l
# order m

# Choose which degrees and order to filter out
minDegree = 1 
maxDegree = 6 
minOrder = 0 
maxOrder = 6

inputfile = 'coeffs.txt' # File containing spherical harmonics created using pyshtools
gridfile = 'expanded_grid.txt' # File name to save the grid values only [To compare or to revert back to SPH later]
filename = 'expanded_grid.vtu' # File name to save the grid in VTU format

print('Reading in SPH data from file "%s"...' % inputfile)
df = pd.read_csv(inputfile, names=['Degree', 'Order', 'Cosine', 'Sine'], dtype=np.float64)

degrees = np.array(df['Degree'])
lastDegree = int(degrees[-1])
coeffs = np.zeros((2, lastDegree + 1, lastDegree + 1))

print('Filtering out degrees outside range: %d-%d \n \t   and orders outside range: %d-%d' % (minDegree, maxDegree, minOrder, maxOrder))
for index, row in df.iterrows():
    if row['Degree'] >= minDegree and row['Degree'] <= maxDegree and row['Order'] >= minOrder and row['Order'] <= maxOrder:
        coeffs[0, int(row['Degree']), int(row['Order'])] = row['Cosine']
        coeffs[1, int(row['Degree']), int(row['Order'])] = row['Sine']

coeffs = pysh.SHCoeffs.from_array(coeffs, normalization='ortho')

nlon = 361
nlat = 181
lons = np.linspace(0, 360, nlon)
lats = np.linspace(-90, 90, nlat)

lats, lons = np.meshgrid(lats, lons)

lons = lons.flatten()
lats = lats.flatten()

outcome = coeffs.expand(lon=lons, lat=lats)

df_out = pd.DataFrame(lons, columns=['Lon'], dtype=np.float64)
df_out['Lat'] = lats
df_out['Outcome'] = outcome

df_out['Lat'] = df_out['Lat'].apply(lambda x: -x)
df_out['Lon'] = df_out['Lon'].apply(lambda x: x - 180)
df_out.sort_values(by=['Lat', 'Lon'], inplace=True)
df.reset_index(drop=True, inplace=True)

lons = np.array(df_out['Lon'], dtype=np.float64)
lats = np.array(df_out['Lat'], dtype=np.float64)
outcome = np.array(df_out['Outcome'], dtype=np.float64)

print('Saving grid values in "%s"...' % gridfile)
np.savetxt(gridfile, outcome, fmt='%.16e', delimiter='\t')

# Connectivity
nnx = nlon
nny = nlat
nelx = nlon - 1 
nely = nlat - 1 
nel = nelx * nely
nnp = nnx * nny 

m = 4 
icon = np.empty((m, nel), dtype=int)
counter = 0 
for ely in range(0, nely):
   for elx in range(0, nelx):
      icon[0, counter] = elx + ely * nnx  
      icon[1, counter] = elx + ely * nnx + 1 
      icon[2, counter] = elx + (ely + 1) * nnx + 1 
      icon[3, counter] = elx + (ely + 1) * nnx 
      counter += 1

print('Saving grid in VTU format in "%s"...' % filename)
vtufile= open(filename, "w")
vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
vtufile.write("<UnstructuredGrid> \n")
vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(nnp,nel))
#####
vtufile.write("<Points> \n")
vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f %10f %10f \n" %(lons[i],lats[i],0.))
vtufile.write("</DataArray>\n")
vtufile.write("</Points> \n")
#####
#####
vtufile.write("<PointData Scalars='scalars'>\n")
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='sph expanded' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %outcome[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("</PointData>\n")

#####
vtufile.write("<Cells>\n")
#--
vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
for iel in range (0,nel):
    vtufile.write("%d %d %d %d\n" %(icon[0,iel],icon[1,iel],icon[2,iel],icon[3,iel]))
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
for iel in range (0,nel):
   vtufile.write("%d \n" %((iel+1)*4))
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
for iel in range (0,nel):
   vtufile.write("%d \n" %9)
vtufile.write("</DataArray>\n")
#--
vtufile.write("</Cells>\n")
#####
vtufile.write("</Piece>\n")
vtufile.write("</UnstructuredGrid>\n")
vtufile.write("</VTKFile>\n")
vtufile.close()

 
