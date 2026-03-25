import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

#we want to calculate the volume transport in the regional boundary files

#bdy file
bdy_path ='/project/6007519/ANHA4-I/BDY/GLORYS2V3/'

bdy_file = bdy_path+'ANHA4_bdy_south_GLORYS2V3_y2002.nc'

bf = xr.open_mfdataset(bdy_file)

#also need the grid file to land mask

grid_file = '/project/6007519/weissgib/ANHA4-I/ANHA4_mesh_mask.nc'

gf = xr.open_mfdataset(grid_file)

e2v = gf['e2v'][0, 0, :].values

#get in the right shape
e2v = np.tile(e2v.reshape(1, 1, 1, 544), (12, 50, 1, 1))

e2v = e2v[:, :, :, 2:544]

gf.close()

#bt = xr.open_mfdataset('/project/6007519/weissgib/ANHA4-I/ANHA4_bathy.nc')
bt = xr.open_mfdataset('/project/6007519/weissgib/ANHA4-I/ANHA4_bathy_etopo1_gebco1_smoothed_coast_corrected_mar10.nc')

bath = bt['Bathymetry'][0, :].values
bath = np.where(bath != 0, 1, 0)
bath = np.tile(bath.reshape(1, 1, 1, 544), (12, 50, 1, 1))

bath = bath[:, :, :, 2:544]

bt.close()

print(bath.shape)

#and the cell thickness, can just use from anha4x

thickness_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/computed_u_thickness.npy'

thickness = np.load(thickness_file)

thickness = thickness[:,0,2:544]

thickness = np.tile(thickness.reshape(1,50,1,542), (12, 1, 1, 1))

print(thickness.shape)

#and now we can calculate the transport

masked_vel = bf['vomecrty']*bath

vol_trans = masked_vel*thickness*e2v

print(vol_trans)

#plt.imshow(vol_trans[0, :, 0, :])
#plt.show()
#exit()

#average the transport over depth and across the line
trans = vol_trans.sum(dim='z')
trans = trans.sum(dim='x')

#convert to sverdrups
trans = trans*0.000001

print(trans)
avg_trans = trans.mean('time_counter')
print(avg_trans.values)

plt.plot(trans[:,0])
#plt.savefig('bdy_glorys_anha4_south.png')
plt.show()

bf.close()



