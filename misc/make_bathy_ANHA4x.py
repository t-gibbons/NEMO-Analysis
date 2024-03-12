"""
Create the ANHA4x bathy and coordinate file

Based off old ANHA4 bathy file, with slice taken from eOrca for extended boundary
"""
import numpy as np
import xarray as xr
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as feature

anha4_file = '/project/6007519/ANHA4-I/ANHA4_bathy_etopo1_gebco1_smoothed_coast_corrected_mar10_Tide.nc'
eorca_file = '/project/6007519/pennelly/eORCA025-I/eORCA025_Bathymetry.nc'

af = xr.open_dataset(anha4_file)
ef = xr.open_dataset(eorca_file)

#get the boundary from old anha4

anha_lon = af['nav_lon'].values
anha_lat = af['nav_lat'].values
anha_bathy = af['Bathymetry'].values
print(anha_bathy.shape)

"""
lon_bd = anha_lon[799, :]
lat_bd = anha_lat[799,:]

print(lon_bd)
print(lat_bd)

#find where that lines up in the eOrca grid
eorca_lon = ef['nav_lon'].values
eorca_lat = ef['nav_lat'].values

lon_start = np.where(eorca_lon == lon_bd[0])
lat_start = np.where(eorca_lat == lat_bd[0])

lon_end = np.where(eorca_lon == lon_bd[-1])
lat_end = np.where(eorca_lat == lat_bd[-1])
"""

#cheating and hardcoding where I think the matching line is in eOrca
#should be j=1009, i=143:686
#want to take a chunk which extends down to j=900

eorca_bathy = ef['Bathymetry'].values
eorca_lon = ef['nav_lon'].values
eorca_lat = ef['nav_lat'].values

eorca_chunk = eorca_bathy[900:1009, 143:687]
lon_chunk = eorca_lon[900:1009, 143:687]
lat_chunk = eorca_lat[900:1009, 143:687]
print(eorca_chunk.shape)

#now need to rotate this chunk so we can properly stitch them together
#use rot90 twice to rotate 180 degrees

r1 = np.rot90(eorca_chunk)
rot_chunk = np.rot90(r1)
print(rot_chunk.shape)

l1 = np.rot90(lon_chunk)
lon_rot = np.rot90(l1)

l2 = np.rot90(lat_chunk)
lat_rot = np.rot90(l2)

#and we need to flip the array
flip_chunk = np.flip(rot_chunk, axis=1)
print(flip_chunk.shape)

flip_lon = np.flip(lon_rot, axis=1)
flip_lat = np.flip(lat_rot, axis=1)

anhax_bathy = np.concatenate((anha_bathy, flip_chunk), axis=0)
print(anhax_bathy.shape)

anhax_lon = np.concatenate((anha_lon, flip_lon), axis=0)
anhax_lat = np.concatenate((anha_lat, flip_lat), axis=0)

#lets try plotting the new anhax on a map
land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)

fig = plt.figure(figsize=(10,9))
ax = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())

ax.add_feature(land_50m, color=[0.8,0.8,0.8])
ax.coastlines(resolution='50m')

p1 = ax.pcolormesh(anhax_lon, anhax_lat, anhax_bathy, transform=ccrs.PlateCarree(), vmin=0, vmax=1000)
ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
cb = plt.colorbar(p1, cax=ax_cb, orientation='vertical')
plt.show()


af.close()
ef.close()
