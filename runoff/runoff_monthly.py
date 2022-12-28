import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

#new runoff
path = '/project/6007519/ANHA4-I/RUNOFF/HydroGFD/'

start_year = 2002
end_year = 2019
data = []

for y in range(start_year, end_year+1):
    files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
    ds = xr.open_mfdataset(files, decode_times=False)
    reference_date = '1/1/'+str(y)
    v = ds['runoff'].values
    ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
    data.append(ds)

new_runoff = xr.concat(data, dim='time_counter')

#old runoff
old_path = '/project/6007519/pmyers/ANHA4-I/RUNOFF/Bamber2012/'
files = []

for y in range(start_year, end_year+1):
    files.append(old_path+'ANHA4_runoff_monthly_combined_Dai_Trenberth_Bamber_y'+str(y)+'.nc')

old_runoff = xr.open_mfdataset(files)

#model grid information for converting units
grid_file = '/project/6007519/weissgib/plotting/ANHA4_mesh_mask.nc'
mesh = xr.open_mfdataset(grid_file)

e1v = mesh['e1v']
e2u = mesh['e2u']
print(e1v)

#and the mask files for getting coastal regions
mask_path = '/project/6007519/weissgib/plotting/runoff_regions_mask.nc'
#mask_path = '/project/6007519/weissgib/plotting/regions_mask.nc'

mask_data = xr.open_mfdataset(mask_path)

#masks = {'ss_mask': 'Siberian Shelf', 'cs_mask': 'Canadian Shelf'}
masks = {'hb_mask': 'Hudson Bay and Greenland', 'caa_mask': 'Canadian Arctic Archipelago', 'bs_mask': 'Bering Strait', 'ec_mask': 'Eastern Coast', 'nc_mask': 'Norway and British Isles'}

#lets make some time series over these regions
for m in masks:
    md = mask_data[m][0,0]
    masked_new_runoff = new_runoff['runoff'].where(md ==2)
    masked_old_runoff = old_runoff['runoff'].where(md ==2)
    #new_timeseries = masked_new_runoff.sum(('x', 'y'))
    #old_timeseries = masked_old_runoff.sum(('x', 'y'))

    e1v_mask = e1v.where(md==2)[0]
    e2u_mask = e2u.where(md==2)[0]

    new_convert = masked_new_runoff*e1v_mask*e2u_mask*0.001
    new_timeseries = new_convert.sum(('x', 'y'))
    old_convert = masked_old_runoff*e1v_mask*e2u_mask*0.001
    old_timeseries = old_convert.sum(('x', 'y'))

    dai_mean = np.mean(old_timeseries.groupby('time_counter.year').mean().values)
    hype_mean = np.mean(new_timeseries.groupby('time_counter.year').mean().values)
    print('Dai annual average mean '+masks[m]+' : '+str(dai_mean))
    print('HYPE annual average mean '+masks[m]+' : '+str(hype_mean))

    nt = new_timeseries.values
    ot = old_timeseries.values

    #lets just repeat the last year of dai to match hype
    last_year = ot[-12:]
    for i in range(0,12): ot = np.concatenate((ot, last_year))

    dn = new_timeseries['time_counter'].values
    n = len(dn)
    ot1 = np.empty(n)
    ot1[:] = np.nan
    for i in range(len(ot)):
        ot1[i] = ot[i]
    
    plt.plot(dn, nt, label='HYPE')
    plt.plot(dn, ot1, label='Dai and Trenberth')
    plt.legend()
    plt.title(masks[m])
    plt.ylabel('runoff (m^3/s)')
    #plt.show()
    plt.savefig('/project/6007519/weissgib/plotting/figs/runoff/'+m+'_runoff_comp_new.png')
    plt.clf()

exit()
#and lets make the same plot but for the whole region
full_new_timeseries = new_runoff['runoff'].sum(('x', 'y'))
full_old_timeseries = old_runoff['runoff'].sum(('x', 'y'))
nt = full_new_timeseries.values
ot = full_old_timeseries.values

last_year = ot[-12:]
for i in range(0,12): ot = np.concatenate((ot, last_year))

dn = new_timeseries['time_counter'].values
n = len(dn)
ot1 = np.empty(n)
ot1[:] = np.nan
for i in range(len(ot)):
    ot1[i] = ot[i]

plt.plot(dn, nt, label='HYPE')
plt.plot(dn, ot1, label='Dai and Trenberth')
#plt.legend()
plt.ylabel('runoff (m^3/s)')

plt.savefig('/project/6007519/weissgib/plotting/figs/runoff/total_runoff_comp_new.png')
plt.clf()

"""
#just take monthly average of the new runoff
#and plot the map for each month

monthly_average = new_runoff.groupby('time_counter.month').mean('time_counter')

mnths = monthly_average['month'].values

for m in range(monthly_average.dims['month']):
    lons = new_runoff['nav_lon'].values
    lats = new_runoff['nav_lat'].values
    mnth = mnths[m]
    values = monthly_average['runoff'][m]
    #north pole stero projection
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=ccrs.NorthPolarStereo()

    fig = plt.figure(figsize=(10, 9))
    ax = plt.subplot(1, 1, 1, projection=projection)
    ax.set_extent([-280, 80, 80, 35], crs=ccrs.PlateCarree())
    #ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='50m')

    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    p1 = ax.pcolormesh(lons, lats, values, transform=ccrs.PlateCarree(), cmap='gist_ncar')
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
    plt.show()
    plt.clf()
"""
new_runoff.close()
old_runoff.close()
mask_data.close()
mesh.close()
