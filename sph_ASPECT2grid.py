import numpy as np
from scipy.special import lpmv
import pandas as pd
import pyshtools as pysh

# degree l
# order m
minDegree = 1 
maxDegree = 50 
minOrder = 0 
maxOrder = 50

inputfile = '../eejit/output-PREM80-5-C-C-A-noVis/geoid_anomaly_SH_coefficients.00001' # SPH file (created by ASPECT) to read 
gridfile = 'prem80_total_1.txt' # File to save the grid values only (For comparison etc.)
filename = 'prem80_total_1.vtu' # File to save the grid in VTU format
rows_to_skip = 2

print('Reading in "%s"...' % inputfile)
df = pd.read_csv(inputfile, skiprows=rows_to_skip, delim_whitespace=True, names=['Degree', 'Order', 'Cosine', 'Sine'], dtype=np.float64)

nlon = 361
nlat = 181
lons = np.linspace(0, 360, nlon)
lats = np.linspace(-90, 90, nlat)

lats, lons = np.meshgrid(lats, lons)

coeffs = np.zeros((2, 41, 41))

for index, row in df.iterrows():
    if row['Degree'] >= minDegree and row['Degree'] <= maxDegree and row['Order'] >= minOrder and row['Order'] <= maxOrder:
        coeffs[0, int(row['Degree']), int(row['Order'])] = row['Cosine']
        coeffs[1, int(row['Degree']), int(row['Order'])] = row['Sine']

coeffs = pysh.SHCoeffs.from_array(coeffs, normalization='ortho')
coeffs.to_file('testcoeffs.txt')
lons = lons.flatten()
lats = lats.flatten()

outcome = coeffs.expand(lon=lons, lat=lats)

df_out = pd.DataFrame(lons, columns=['Lon'], dtype=np.float64)
df_out['Lat'] = lats
df_out['Outcome'] = outcome

df_out['Lon'] = df_out['Lon'].apply(lambda x: x - 180)
df_out.sort_values(by=['Lat', 'Lon'], inplace=True)
df.reset_index(drop=True, inplace=True)

lons = np.array(df_out['Lon'], dtype=np.float64)
lats = np.array(df_out['Lat'], dtype=np.float64)
outcome = np.array(df_out['Outcome'], dtype=np.float64)

print('Saving values on grid in "%s"...' % gridfile)
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

 
