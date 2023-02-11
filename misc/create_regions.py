"""
create_regions.py
author: Tahya Weiss-Gibbons (weissgib@ualberta.ca)
date: June 2021

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
#mesh_file = '/mnt/storage4/tahya/model_files/ANHA12_mesh_zgr.nc'
mesh_data = xr.open_dataset(mesh_file)

x = mesh_data.dims['x']
y = mesh_data.dims['y']
z = mesh_data.dims['z']

nav_lon = mesh_data['nav_lon']
nav_lat = mesh_data['nav_lat']
nav_lev = mesh_data['mbathy']

#want to define the lat/lon boundaries of the regions we want

#old definition
#can_basin = {'lon':[-161, -120], 'lat': [70, 80]}
can_basin = {'lon':[-170, -130], 'lat':[70.5,80.5]}
eur_basin = {'lon':[80, 168], 'lat':[78, 88]}
caa = {'lon':[-127, -80], 'lat':[66, 80]}
hudson_bay = {'lon':[-95, -66], 'lat':[50, 64]}
can_shelf = {'lon':[-171, -127], 'lat':[68, 76]}
siberian_shelf = {'lon': [132,177], 'lat':[68,76]}
bering_strait = {'lon':[-180,-158], 'lat': [63,71]}
cntrl_arctic = {'lon':[-122,22], 'lat':[84,89]}
davis_strait = {'lon': [-78,-53], 'lat':[66,74]}
lab_sea = {'lon':[-51,-23], 'lat':[43,59]}
nares_strait = {'lon': [-80, -34], 'lat': [74, 84]}
fram_strait = {'lon': [-41, -4], 'lat': [62, 81]}
lab_current = {'lon': [-70, -48], 'lat': [41, 60]}
spna = {'lon': [0,60], 'lat': [47, 65]}

#now go through the grid and make surface masks

mesh_data['cb_mask'] = mesh_data['tmask'].copy()
mesh_data['eb_mask'] = mesh_data['tmask'].copy()
mesh_data['caa_mask'] = mesh_data['tmask'].copy()
mesh_data['hb_mask'] = mesh_data['tmask'].copy()
mesh_data['ss_mask'] = mesh_data['tmask'].copy()
mesh_data['cs_mask'] = mesh_data['tmask'].copy()
mesh_data['bs_mask'] = mesh_data['tmask'].copy()
mesh_data['ca_mask'] = mesh_data['tmask'].copy()
mesh_data['ds_mask'] = mesh_data['tmask'].copy()
mesh_data['ls_mask'] = mesh_data['tmask'].copy()
mesh_data['ns_mask'] = mesh_data['tmask'].copy()
mesh_data['fs_mask'] = mesh_data['tmask'].copy()
mesh_data['lc_mask'] = mesh_data['tmask'].copy()
mesh_data['spna_mask'] = mesh_data['tmask'].copy()

mesh_data['all_masks'] = mesh_data['tmask'].copy()

for i in range(x):
    for j in range(y):
        ln = nav_lon[j,i]
        lt = nav_lat[j,i]
        lv = nav_lev[0,j,i] 
		
        if can_basin['lon'][0] <= ln <= can_basin['lon'][1]:
            if can_basin['lat'][0] <= lt <= can_basin['lat'][1]:
                if lv >= 28:
                    mesh_data['cb_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 2

        if eur_basin['lon'][0] <= ln <= eur_basin['lon'][1]:
            if eur_basin['lat'][0] <= lt <= eur_basin['lat'][1]:
                if lv >= 35:
                    mesh_data['eb_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 3

        if caa['lon'][0] <= ln <= caa['lon'][1]:
            if caa['lat'][0] <= lt <= caa['lat'][1]:
                if lv <= 31:
                    mesh_data['caa_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 4

        if hudson_bay['lon'][0] <= ln <= hudson_bay['lon'][1]:
            if hudson_bay['lat'][0] <= lt <= hudson_bay['lat'][1]:
                mesh_data['hb_mask'][:,:,j,i] = 2
                mesh_data['all_masks'][:,:,j,i] = 5

        if siberian_shelf['lon'][0] <= ln <= siberian_shelf['lon'][1]:
            if siberian_shelf['lat'][0] <= lt <= siberian_shelf['lat'][1]:
                if lv <= 20:
                    mesh_data['ss_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 6

        if can_shelf['lon'][0] <= ln <= can_shelf['lon'][1]:
            if can_shelf['lat'][0] <= lt <= can_shelf['lat'][1]:
                if lv <= 20:
                    mesh_data['cs_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 7

        if bering_strait['lon'][0] <= ln <= bering_strait['lon'][1]:
            if bering_strait['lat'][0] <= lt <= bering_strait['lat'][1]:
                mesh_data['bs_mask'][:,:,j,i] = 2
                mesh_data['all_masks'][:,:,j,i] = 8

        if cntrl_arctic['lon'][0] <= ln <= cntrl_arctic['lon'][1]:
            if cntrl_arctic['lat'][0] <= lt <= cntrl_arctic['lat'][1]:
                mesh_data['ca_mask'][:,:,j,i] = 2
                mesh_data['all_masks'][:,:,j,i] = 9

        if davis_strait['lon'][0] <= ln <= davis_strait['lon'][1]:
            if davis_strait['lat'][0] <= lt <= davis_strait['lat'][1]:
                if lv >= 22:
                    mesh_data['ds_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 10

        if lab_sea['lon'][0] <= ln <= lab_sea['lon'][1]:
            if lab_sea['lat'][0] <= lt <= lab_sea['lat'][1]:
                if lv >= 35:
                    mesh_data['ls_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 11

        if lab_current['lon'][0] <= ln <= lab_current['lon'][1]:
            if lab_current['lat'][0] <= lt <= lab_current['lat'][1]:
                if lv <= 35:
                    mesh_data['lc_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 12

        if nares_strait['lon'][0] <= ln <= nares_strait['lon'][1]:
            if nares_strait['lat'][0] <= lt <= nares_strait['lat'][1]:
                mesh_data['ns_mask'][:,:,j,i] = 2
                mesh_data['all_masks'][:,:,j,i] = 13

        if fram_strait['lon'][0] <= ln <= fram_strait['lon'][1]:
            if fram_strait['lat'][0] <= lt <= fram_strait['lat'][1]:
                if lv <= 45:
                    mesh_data['fs_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 14

        if spna['lon'][0] <= ln <= spna['lon'][1]:
            if spna['lat'][0] <= lt <= spna['lat'][1]:
                if lv >= 35:
                    mesh_data['spna_mask'][:,:,j,i] = 2
                    mesh_data['all_masks'][:,:,j,i] = 15

        

#remask the data with the bathymetry now, since we ignored that when finding our bounding boxes
mesh_data = mesh_data.where(mesh_data.tmask ==1)

#lets write the data we have to a netcdf file

mesh_data.to_netcdf('new_region_mask.nc')

