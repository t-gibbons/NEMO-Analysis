import glob
import numpy as np
import xarray as xr
import netCDF4 as nc

#constants

rho = 1025 #kg/m^3
c_p = 3850 #K/(kgC)
runid = 'EPM161'

#paths
path = "/project/6007519/weissgib/ANHA4/ANHA4-"+runid+"-S/"
grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'
output_path = '/project/6007519/weissgib/plotting/heat/'

#first read the model files
mdl_files = glob.glob(path+'ANHA4-'+runid+'*_gridT.nc')

#remove the last file from the list
#for a lot of the model runs, the last file is blank so just going to remove it
#last_file = path+'ANHA4-'+runid+'_y2018m11d11_gridT.nc' #just set last file for now, want to fix this
#mdl_files.remove(last_file)

d = xr.open_mfdataset(mdl_files, concat_dim='time_counter', data_vars='minimal', coords='minimal', compat='override')

mesh = nc.Dataset(grid_file)

mask = np.array(mesh.variables['tmask'])

mesh.close()

#mask the data

d.coords['mask'] = (('deptht', 'y_grid_T', 'x_grid_T'), mask[0,:,:,:])
d = d.where(d.mask == 1)

d = d.where(d['deptht'] < 100.0, drop=True) #depth to calculate the heat content over

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

tmp = rho*c_p*(d-1.8)

x = d.sizes['x_grid_T']
y = d.sizes['y_grid_T']

#get dz on the grid
dz_grid = np.tile(dz[:,np.newaxis,np.newaxis], (1,y,x))

t1 = tmp*dz_grid

heat_content = t1.sum(dim='deptht', skipna=True)

heat_content.to_netcdf(output_path+runid+'_heat_content.nc')
