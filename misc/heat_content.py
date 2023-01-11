import glob
import numpy as np
import xarray as xr
import netCDF4 as nc

#constants

rho = 1025 #kg/m^3
c_p = 3850 #K/(kgC)

#paths
#path = "/project/6000276/weissgib/ANHA4/ANHA4-ETW101-S/"
path = "/project/6007519/pmyers/ANHA4/ANHA4-EPM101-S/"
grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'

#first read the model files
mdl_files = glob.glob(path+'ANHA4-EPM101*_gridT.nc')

#remove the last file from the list
last_file = path+'ANHA4-EPM101_y2019m04d10_gridT.nc' #want to figure this out automatically ideally
mdl_files.remove(last_file)

d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

mesh = nc.Dataset(grid_file)

mask = np.array(mesh.variables['tmask'])

mesh.close()

#mask the data

d.coords['mask'] = (('deptht', 'y_grid_T', 'x_grid_T'), mask[0,:,:,:])
d = d.where(d.mask == 1)

d = d.where(d['deptht'] < 100.0, drop=True)

d = d['votemper']

full_depth = list(d['deptht'].values)

#calculate the weights of the depth levels
n = len(d['deptht'])
weight = np.zeros(n)
dz = np.zeros(n)
dd = d['deptht'][n-1]
for i in range(n):
    if i == 0:
        weight[i] = d['deptht'][i]/dd
        dz[i] = d['deptht'][i]
    else:
        weight[i] = (d['deptht'][i] - d['deptht'][i-1])/dd
        dz[i] = d['deptht'][i] - d['deptht'][i-1]

tmp = (-1.8-d)/-1.8

x = d.sizes['x_grid_T']
y = d.sizes['y_grid_T']

#get dz on the grid
dz_grid = np.tile(dz[:,np.newaxis,np.newaxis], (1,y,x))

t1 = tmp*dz_grid

heat_content = t1.sum(dim='deptht', skipna=True)
print(heat_content)

heat_content.to_netcdf('/project/6000276/weissgib/plotting/EPM101_heat_content.nc')
