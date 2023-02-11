"""
create_regions.py
author: Tahya Weiss-Gibbons (weissgib@ualberta.ca)
date: June 2021
updated: January 2022

creates mask files for different regions on the ANHA4 grid
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

#first lets read in the mesh grid file
mesh_file = '/mnt/storage4/tahya/model_files/ANHA4_mesh_mask.nc'
mesh_data = xr.open_dataset(mesh_file)

x = mesh_data.dims['x']
y = mesh_data.dims['y']
z = mesh_data.dims['z']

nav_lon = mesh_data['nav_lon']
nav_lat = mesh_data['nav_lat']
nav_lev = mesh_data['mbathy']

#want to define the lat/lon boundaries of the regions we want

hb_coast = {'lon': [-83.77, -15.52], 'lat': [58.49, 84]}
hudson_bay = {'lon':[-98, -76], 'lat':[50, 66]}
caa = {'lon':[-127, -80], 'lat':[66, 80]}
#bs_coast = {'lon':[248.91,124.27], 'lat':[64.7,80.5]}
bs_coast = {'lon':[-167, -125], 'lat': [59, 72]}
#east_coast = {'lon': [31.43, 131.81], 'lat':[64.35, 76.76]}
norway_coast = {'lon': [-18.14, 41], 'lat':[52.25,71.66]}
bs_east = {'lon':[140, 190], 'lat': [66,77]}
laptev_sea = {'lon': [104, 142], 'lat': [68, 80]}
kara_sea = {'lon': [53, 104], 'lat': [65, 79]}

#now go through the grid and make surface masks

mesh_data['hb_mask'] = mesh_data['tmask'].copy()
mesh_data['caa_mask'] = mesh_data['tmask'].copy()
mesh_data['bs_mask'] = mesh_data['tmask'].copy()
mesh_data['ec_mask'] = mesh_data['tmask'].copy()
#mesh_data['nc_mask'] = mesh_data['tmask'].copy()
mesh_data['bs_east_mask'] = mesh_data['tmask'].copy()
mesh_data['laptev_mask'] = mesh_data['tmask'].copy()
mesh_data['kara_mask'] = mesh_data['tmask'].copy()
mesh_data['all_masks'] = mesh_data['tmask'].copy()

for i in range(x):
    for j in range(y):
        ln = nav_lon[j,i]
        lt = nav_lat[j,i]
        lv = nav_lev[0,j,i] 
	
        """	
        if hb_coast['lon'][0] <= ln <= hb_coast['lon'][1]:
            if hb_coast['lat'][0] <= lt <= hb_coast['lat'][1]:
                if lv <= 28:
                    mesh_data['hb_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 2
        """
        if hudson_bay['lon'][0] <= ln <= hudson_bay['lon'][1]:
            if hudson_bay['lat'][0] <= lt <= hudson_bay['lat'][1]:
                mesh_data['hb_mask'][:,:,j,i] = 2
                mesh_data['all_masks'][:,:,j,i] = 2

        if caa['lon'][0] <= ln <= caa['lon'][1]:
            if caa['lat'][0] <= lt <= caa['lat'][1]:
                if lv <= 25:
                    mesh_data['caa_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 3

        #bering strait section is a special case
        #it crosses the lon boundary in -180 to 180 case
        #so going to convert to 0 360 to avoid the problem

        #dln = ((360 + (ln % 360)) % 360)
        if bs_coast['lon'][0] <= ln <= bs_coast['lon'][1]: 
            if bs_coast['lat'][0] <= lt <= bs_coast['lat'][1]:
               if lv <= 25:
                   mesh_data['bs_mask'][:,:,j,i] = 2
                   mesh_data['all_masks'][:,:,j,i] = 4

        if bs_east['lon'][0] <= ln <= bs_east['lon'][1]:
            if bs_east['lat'][0] <= lt <= bs_east['lat'][1]:
                if lv <= 25:
                    mesh_data['bs_east_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 5

        if laptev_sea['lon'][0] <= ln <= laptev_sea['lon'][1]:
            if laptev_sea['lat'][0] <= lt <= laptev_sea['lat'][1]:
                if lv <= 25:
                    mesh_data['laptev_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 6

        if kara_sea['lon'][0] <= ln <= kara_sea['lon'][1]:
            if kara_sea['lat'][0] <= lt <= kara_sea['lat'][1]:
                if lv <= 25:
                    mesh_data['kara_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 7


#        if norway_coast['lon'][0] <= ln <= norway_coast['lon'][1]:
#            if norway_coast['lat'][0] <= lt <= norway_coast['lat'][1]:
#                if lv <= 25:
#                    mesh_data['nc_mask'][:,:,j,i] = 2
#                    mesh_data['all_masks'][:,:,j,i] = 8



#remask the data with the bathymetry now, since we ignored that when finding our bounding boxes
mesh_data = mesh_data.where(mesh_data.tmask ==1)

#lets write the data we have to a netcdf file

mesh_data.to_netcdf('runoff_temp_regions_mask_all_mask.nc')

