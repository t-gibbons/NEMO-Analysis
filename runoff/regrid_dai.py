import numpy as np
import pandas as pd
import xarray as xr
import geopy.distance
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

#get the model data for the new grid

model_path = '/mnt/storage6/clark/eORCA025/eORCA025_mesh_mask.nc'

model_data = xr.open_mfdataset(model_path)
model_lat = model_data['nav_lat'].values
model_lon = model_data['nav_lon'].values
e1 = model_data['e1t'][0].values
e2 = model_data['e2t'][0].values

model_data.close()

#also need the runoff mask
mask_path = '/mnt/storage5/clark/eORCA025/Runoff/eORCA025_ReNat_HydroGFD_HBC_runoff_monthly_y1996.nc'

mask_data = xr.open_mfdataset(mask_path)
runoff_mask = mask_data['socoefr'].values

print(runoff_mask.shape)

mask_data.close()

dai_path = '/mnt/storage6/tahya/data/dai_2021_runoff/'

dai = xr.open_mfdataset(dai_path+'coastal-stns-Vol-monthly.updated-May2019.nc', decode_times=False)

print(dai['time'].values)

dai['time'] = pd.to_datetime(dai['time'].values, format="%Y%m")
print(dai['time'])

#now seperate output files into years
years = ['2000', '2018']
for year in years:
    times = pd.date_range(start='1/1/'+year, periods=12, freq='1M')
    print(times)

    #construct the regridded dai dataset
    runoff = np.zeros((12, model_data.sizes['y'], model_data.sizes['x']))
    rg_dai = xr.Dataset(
        data_vars=dict(
            runoff=(["t","y", "x"], runoff),
            socoefr=(["y","x"], runoff_mask),
            time_counter=(["t"], times),
            nav_lon=(["y", "x"], model_lon),
            nav_lat=(["y", "x"], model_lat),
        ),
        attrs=dict(description="Dai 2021 Runoff Data: Regridded eOrca grid"),
    )

    #find the closest model grid cell to each station
    for s in dai['station'].values:
        ln = dai['lon'].sel(station=s).values
        lt = dai['lat'].sel(station=s).values
        abslat = np.abs(model_lat-lt)
        abslon = np.abs(model_lon-ln)

        c = np.maximum(abslon,abslat)
        c_mask = np.ma.masked_where(runoff_mask == 0, c)    #mask the array so we only consider grid cells on the coast
        y, x = np.where(c_mask == np.min(c_mask))
    
        #sometimes there are two closest grid cells
        #just take the first one in this case
        if y.shape != (1,):
            y = y[:-1]
            x = x[:-1]
        mdl_ln = model_lon[y[0],x[0]]
        mdl_lt = model_lat[y[0],x[0]]

        flow = dai['FLOW'].sel(station=s,time=dai.time.dt.year.isin([int(year)])).values
        grid_flow = rg_dai['runoff'][:,y,x].values
        total = np.add(flow, grid_flow[:,0,0])
        total = total[..., np.newaxis, np.newaxis]

        #convert the units
        #going from m^3/s to kg/s/m^2
        #density of water 1028 kg/m^3
        convert_total = (total*1028)/(e1[y,x]*e2[y,x])
        convert_total[np.isnan(convert_total)] = 0

        rg_dai['runoff'].loc[:,y,x] += convert_total

    """
    #test if locations approimately match
    projection = ccrs.Mercator()
    fig = plt.figure(figsize=(10,9))
    ax = plt.subplot(1, 1, 1, projection=projection)
    ax.coastlines(resolution='50m')

    ax.plot(rg_dai['nav_lon'], rg_dai['nav_lat'], rg_dai['runoff'][0], transform=ccrs.PlateCarree())
    ax.scatter(dai['lon'], dai['lat'], transform=ccrs.PlateCarree(), s=50)
    plt.show()
    #plt.savefig('dai_eorca_'+year+'_stat_loc.png')
    """

    print(rg_dai)
    rg_dai.to_netcdf('dai_remapped/dai_eorca_'+year+'.nc')

model_data.close()
dai.close()
