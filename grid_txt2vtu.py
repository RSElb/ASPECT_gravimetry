import numpy as np
import pandas as pd
import sys as sys

try:
    inputfile = str(sys.argv[1])
    filename = outputfile = str(sys.argv[2])
    dataname = str(sys.argv[3])
except:
    print('Please enter one input (.txt file) and one output (.vtu) file and the name of the data (and optionally a scaling factor)!')

try:
    scaling_factor = float(sys.argv[4])
except:
    scaling_factor = 1.

print('Reading in grid (.txt) data from "%s"...' % inputfile)
df = pd.read_csv(inputfile, delim_whitespace=True, names=['Value'], dtype=np.float64)

outcome = np.array(df['Value'], dtype=np.float64)

nlon = 361
nlat = 181
lons = np.linspace(-180, 180, nlon)
lats = np.linspace(-90, 90, nlat)

lats, lons = np.meshgrid(lats, lons)
lons = lons.flatten()
lats = lats.flatten()

df_out = pd.DataFrame(lons, columns=['Lon'], dtype=np.float64)
df_out['Lat'] = lats

df_out.sort_values(by=['Lat', 'Lon'], inplace=True)
df.reset_index(drop=True, inplace=True)

lons = np.array(df_out['Lon'], dtype=np.float64)
lats = np.array(df_out['Lat'], dtype=np.float64)

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
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='%s' Format='ascii'> \n" % dataname)
for i in range(0,nnp):
    vtufile.write("%10f \n" %(outcome[i]/scaling_factor))
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
