import numpy as np
import pandas as pd
import sys as sys
from scipy.interpolate import griddata
import math as math

def lon_correction(lon):
   if lon > 180:
	   return lon - 360
   else:
      return lon

def lat_correction(lat):
	return 90 - lat 

try:
   inputfile = sys.argv[1]
except:
   inputfile = 'gravity-00000'

try:
   scaling_factor = float(sys.argv[2])
except:
   scaling_factor = 1.

try:
   remove = float(sys.argv[3])
except:
   remove = 0.

try:
   nlon = int(sys.argv[4])
   nlat = int(sys.argv[5])
except:
   nlon = 361
   nlat = 181

Re = 6371e3 # Radius Earth
Rs = 6596e3 # Radius satellite
Rs = Re # geoid is measured at the surface...

lons = np.linspace(-180, 180, nlon)
lats = np.linspace(-90, 90, nlat)

print('Reading in file "%s" and saving in VTU format in file "gravity_map.vtu" and "gravity_sphere.vtu"' % (inputfile))

df = pd.read_csv(inputfile, skiprows=30, delim_whitespace=True, dtype=np.float64, names=['pos_sat_r', 'pos_sat_phi', 'pos_sat_theta', 'pos_sat_x', 'pos_sat_y', 'pos_sat_z', 'grav_x', 'grav_y', 'grav_z', 'grav_norm', 'grav_theory', 'grav_pot', 'grav_pot_theory', 'grav_anom_x', 'grav_anom_y', 'grav_anom_z', 'grav_anom_norm', 'grav_grad_xx', 'grav_grad_yy', 'grav_grad_zz', 'grav_grad_xy', 'grav_grad_xz', 'grav_grad_yz', 'grav_grad_xx_theory', 'grav_grad_yy_theory', 'grav_grad_zz_theory', 'grav_grad_xy_theory', 'grav_grad_xz_theory', 'grav_grad_yz_theory'])

r = 6371.e3
df['x'] = r * np.sin(np.radians(df['pos_sat_theta'])) * np.cos(np.radians(df['pos_sat_phi']))
df['y'] = r * np.sin(np.radians(df['pos_sat_theta'])) * np.sin(np.radians(df['pos_sat_phi']))
df['z'] = r * np.cos(np.radians(df['pos_sat_theta']))

df['pos_sat_phi'] = df['pos_sat_phi'].apply(lambda x: lon_correction(x))
df['pos_sat_theta'] = df['pos_sat_theta'].apply(lambda x: lat_correction(x))

df.sort_values(by=['pos_sat_theta', 'pos_sat_phi'], inplace=True)
df.reset_index(drop=True, inplace=True)

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

x = np.array(df['x'], dtype=np.float32)
y = np.array(df['y'], dtype=np.float32)
z = np.array(df['z'], dtype=np.float32)
lons = np.array(df['pos_sat_phi'], dtype=np.float32)
lats = np.array(df['pos_sat_theta'], dtype=np.float32)
sat_x = np.array(df['pos_sat_x'], dtype=np.float32)
sat_y = np.array(df['pos_sat_y'], dtype=np.float32)
sat_z = np.array(df['pos_sat_x'], dtype=np.float32)
grav_x = np.array(df['grav_x'], dtype=np.float32)
grav_y = np.array(df['grav_y'], dtype=np.float32)
grav_z = np.array(df['grav_z'], dtype=np.float32)
grav_norm = np.array(df['grav_norm'], dtype=np.float64)
grav_theory = np.array(df['grav_theory'], dtype=np.float32)
grav_pot = np.array(df['grav_pot'], dtype=np.float64)
grav_pot_theory = np.array(df['grav_pot_theory'], dtype=np.float32)
grav_anom_x = np.array(df['grav_anom_x'], dtype=np.float32)
grav_anom_y = np.array(df['grav_anom_y'], dtype=np.float32)
grav_anom_z = np.array(df['grav_anom_z'], dtype=np.float32)
grav_anom_norm = np.array(df['grav_anom_norm'], dtype=np.float32)
grav_grad_xx = np.array(df['grav_grad_xx'], dtype=np.float32)
grav_grad_yy = np.array(df['grav_grad_yy'], dtype=np.float32)
grav_grad_zz = np.array(df['grav_grad_zz'], dtype=np.float32)
grav_grad_xy = np.array(df['grav_grad_xy'], dtype=np.float32)
grav_grad_xz = np.array(df['grav_grad_xz'], dtype=np.float32)
grav_grad_yz = np.array(df['grav_grad_yz'], dtype=np.float32)
grav_grad_xx_theory = np.array(df['grav_grad_xx_theory'], dtype=np.float32)
grav_grad_yy_theory = np.array(df['grav_grad_yy_theory'], dtype=np.float32)
grav_grad_zz_theory = np.array(df['grav_grad_zz_theory'], dtype=np.float32)
grav_grad_xy_theory = np.array(df['grav_grad_xy_theory'], dtype=np.float32)
grav_grad_xz_theory = np.array(df['grav_grad_xz_theory'], dtype=np.float32)
grav_grad_yz_theory = np.array(df['grav_grad_yz_theory'], dtype=np.float32)

mean = np.mean(grav_norm)
grav_anomaly = grav_norm - mean

grav_pot_surface = grav_pot * (Rs / Re)
mean_for_geoid = np.mean(grav_pot_surface)

gravity_for_grid = (grav_norm / scaling_factor) - remove
geoid_for_grid = (grav_pot_surface - mean_for_geoid) / -9.81

print('Saving gravity grid in "gravity_grid.txt"...')
print(type(geoid_for_grid[0]))
np.savetxt('gravity_grid.txt', gravity_for_grid, fmt='%.32e')

print('Saving geoid grid in "geoid_grid.txt"...')
np.savetxt('geoid_grid.txt', geoid_for_grid, fmt='%.32e')

filename = 'gravity_map.vtu'
print('Saving data in VTU format in "%s"...' % filename)
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
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_x' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_x[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_y' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_y[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_z' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_z[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_norm' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %(grav_norm[i]/scaling_factor - remove))
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_potential' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_pot[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_potential_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_pot_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_x' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_x[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_y' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_y[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_z' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_z[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_norm' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_norm[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xx' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xx[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yy' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yy[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_zz' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_zz[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xy' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xy[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xz' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xz[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yz' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yz[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xx_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xx_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yy_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yy_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_zz_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_zz_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xy_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xy_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xz_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xz_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yz_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yz_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_calculated' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anomaly[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='geoid' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" % ((grav_pot_surface[i] - mean_for_geoid) / -9.81))
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

# Connectivity
nnx = nlon
nny = nlat
nelx = nlon 
nely = nlat - 1
nel = nelx * nely
nnp = nnx * nny

m = 4
icon = np.empty((m, nel), dtype=int)
counter = 0
for ely in range(0, nely):
   for elx in range(0, nelx):
      icon[0, counter] = elx + ely * nnx	
      if elx == nlon-1:
         icon[1, counter] = elx + (ely - 1) * nnx + 1
         icon[2, counter] = elx + (ely) * nnx + 1
      else:
         icon[1, counter] = elx + ely * nnx + 1
         icon[2, counter] = elx + (ely + 1) * nnx + 1
      icon[3, counter] = elx + (ely + 1) * nnx
      counter += 1

filename = 'gravity_sphere.vtu'
vtufile= open(filename, "w")
vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
vtufile.write("<UnstructuredGrid> \n")
vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(nnp,nel))
#####
vtufile.write("<Points> \n")
vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f %10f %10f \n" %(x[i],y[i],z[i]))
vtufile.write("</DataArray>\n")
vtufile.write("</Points> \n")
#####
#####
vtufile.write("<PointData Scalars='scalars'>\n")
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_x' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_x[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_y' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_y[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_z' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_z[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_norm' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %(grav_norm[i]/scaling_factor - remove))
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_potential' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_pot[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_x' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_x[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_y' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_y[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_z' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_z[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_anomaly_norm' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_anom_norm[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xx' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xx[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yy' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yy[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_zz' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_zz[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xy' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xy[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xz' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xz[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yz' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yz[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xx_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xx_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yy_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yy_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_zz_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_zz_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xy_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xy_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_xz_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_xz_theory[i])
vtufile.write("</DataArray>\n")
#--
vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='gravity_gradient_yz_theory' Format='ascii'> \n")
for i in range(0,nnp):
   vtufile.write("%10f \n" %grav_grad_yz_theory[i])
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
