import numpy as np
import xarray as xr
import pandas as pd
import scipy.io as sio
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.feature as feature

def read_ardat():
    #lets read in the spreadsheet data

    s_file = '/mnt/storage6/tahya/obs/ARDAT.mat'

    ss = sio.loadmat(s_file, squeeze_me=True)

    #need the lat, lon, temp
    lat = ss['lat']
    lon = ss['long']
    temp = ss['temperature']

    return lat,lon,temp

def region_timeseries():
    
    lat, lon, temp = read_ardat()

    #how many locations do we actually have data for?
    temp1 = temp[~np.all(temp == 0, axis=2)]
    print(temp1.shape)

    #region we want to compare with
    region =  {'lon':[-167, -125], 'lat': [59, 72]}

    coord = lat.shape
    print(coord)

    river_data = []
    for i in range(coord[0]):
        for j in range(coord[1]):
           if region['lon'][0] <= lon[i,j] <= region['lon'][1]:
               if region['lat'][0] <= lat[i,j] <= region['lat'][1]:
                   river_data.append(temp[i,j,:])

    #now should have a list of arrays length 12
    #want to average them for each month

    river_data = np.asarray(river_data)
    print(river_data.shape)

    #most of these locations have no data
    #lets remove the values that are all zero

    river_data = river_data[~np.all(river_data == 0, axis=1)]
    print(river_data.shape)

    mean_obs = river_data.mean(axis=0)

    #now we want the same information from the A-HYPE data to compare
    #new runoff
    path = '/mnt/storage4/tahya/runoff/runoff_temp_files/'

    start_year = 2002
    end_year = 2018
    data = []

    for y in range(start_year, end_year+1):
        if y == 2011: continue
        files = path+'ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y'+str(y)+'.nc'
        ds = xr.open_mfdataset(files, decode_times=False)
        reference_date = '1/1/'+str(y)
        v = ds['rotemper'].values
        ds['time_counter'] = pd.date_range(start=reference_date, periods=ds.sizes['time_counter'], freq='MS')
        data.append(ds)

    new_runoff = xr.concat(data, dim='time_counter')
    new_runoff = new_runoff.where(new_runoff['rotemper'] > -999)
    new_runoff = new_runoff.where(new_runoff['runoff'] != 0)

    #and lets use the bs_east mask, same definition as used for the observations
    mask_path = '/mnt/storage4/tahya/runoff/runoff_temp_regions_mask.nc'

    mask_data = xr.open_mfdataset(mask_path)

    md = mask_data['bs_mask'][0,0]
    masked_runoff = new_runoff['rotemper'].where(md ==2)
    timeseries = masked_runoff.mean(('x', 'y'))

    hype_mean = timeseries.groupby('time_counter.month').mean()

    print(hype_mean)

    nt = hype_mean.values

    dn = hype_mean['month'].values

    plt.plot(dn,nt,label='A-HYPE')
    plt.plot(dn,mean_obs,label='ARDAT')
    plt.legend()
    plt.title('Monthly Average River Water Temperatures for Mackenzie Region')
    plt.ylabel('Temperature (C)')
    plt.xlabel('Month')

    plt.savefig('model_obs_mackenzie_water_temp.png')

def plot_ardat():
    lat, lon, temp = read_ardat()

    #just need to average the temperature data across the year
    temp = temp.mean(axis=2)
    print(temp.shape)

    #and convert all the zero points to nan
    temp[temp == 0] = np.nan

    #map projection
    land_50m = feature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5)
    projection=ccrs.NorthPolarStereo()

    fig = plt.figure(figsize=(10, 9))
    ax = plt.subplot(1, 1, 1, projection=projection)

    ax.set_extent([-180,180,90,60], crs=ccrs.PlateCarree())
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
    
    p1 = ax.scatter(lon, lat, c=temp, transform=ccrs.PlateCarree(), cmap='Reds')
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel('River water temp')
    ax.gridlines()
    #plt.show()
    plt.savefig('obs_temp_locations.png')    

if __name__ == "__main__":
    region_timeseries()

