"""
Create the ANHA4x bathy and coordinate file

Based off old ANHA4 bathy file, with slice taken from eOrca for extended boundary

author: Tahya Weiss-Gibbons (weissgib@ualberta.ca)
"""
import re
import glob
import pathlib
import numpy as np
import xarray as xr
import netCDF4 as nc


def make_bathy():
    anha4_file = '/project/6007519/ANHA4-I/ANHA4_bathy_etopo1_gebco1_smoothed_coast_corrected_mar10_Tide.nc'
    eorca_file = '/project/6007519/pennelly/eORCA025-I/eORCA025_Bathymetry.nc'

    af = xr.open_dataset(anha4_file)
    ef = xr.open_dataset(eorca_file)

    #get the boundary from old anha4

    anha_lon = af['nav_lon'].values
    anha_lat = af['nav_lat'].values
    anha_bathy = af['Bathymetry'].values
    print(anha_bathy.shape)

    eorca_bathy = ef['Bathymetry'].values
    eorca_lon = ef['nav_lon'].values
    eorca_lat = ef['nav_lat'].values

    """
    #this is how figured out where to cut eorca to match anha
    lon_bd = anha_lon[799, :]
    lat_bd = anha_lat[799,:]

    #find where that lines up in the eOrca grid

    lon_start = np.where(eorca_lon == lon_bd[0])
    lat_start = np.where(eorca_lat == lat_bd[0])

    lon_end = np.where(eorca_lon == lon_bd[-1])
    lat_end = np.where(eorca_lat == lat_bd[-1])
    """

    #cheating and hardcoding where I think the matching line is in eOrca
    #should be j=1009, i=143:686, found using output above
    #want to take a chunk which extends down to j=898

    eorca_bathy = ef['Bathymetry'].values
    eorca_lon = ef['nav_lon'].values
    eorca_lat = ef['nav_lat'].values

    eorca_chunk = eorca_bathy[898:1009, 143:687]
    lon_chunk = eorca_lon[898:1009, 143:687]
    lat_chunk = eorca_lat[898:1009, 143:687]
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

    anhax_bathy = np.concatenate((anha_bathy, rot_chunk), axis=0)
    print(anhax_bathy.shape)

    anhax_lon = np.concatenate((anha_lon, lon_rot), axis=0)
    anhax_lat = np.concatenate((anha_lat, lat_rot), axis=0)

    af.close()
    ef.close()

def make_coords():
    anha4_file = '/project/6007519/ANHA4-I/ANHA4_coordinates.nc'
    eorca_file = '/project/6007519/pennelly/eORCA025-I/eORCA025_coordinates.nc'

    af = xr.open_dataset(anha4_file)
    ef = xr.open_dataset(eorca_file)

    #we know the dims will be y=911 and x=544 for the new anhax
    anhax = xr.Dataset(None, coords=dict(y=range(911), x=range(544)))

    for v in list(af.keys()):
        print(v)
        anha_var = af[v].values
        ev = ef[v].values

        #get the variable chunk to add from eorca
        eorca_chunk = ev[898:1009, 143:687]

        #now rotate to match the anha grid

        r1 = np.rot90(eorca_chunk)
        rot_chunk = np.rot90(r1)

        #and concatenate with anha
        anhax_var = np.concatenate((anha_var, rot_chunk), axis=0)
        print(anhax_var.shape)

        anhax[v] = (('y','x'), anhax_var)

    print(anhax)
    
    anhax.to_netcdf('ANHAx_coordinates.nc')

    af.close()
    ef.close()

if __name__ == "__main__":
   make_bathy()

