import glob
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

def max_timeseries():
    path = '/project/6007519/ANHA4-I/ATMFORCING/ERA5/Ready/'

    atmo_files = glob.glob(path+'NEMO_ERA5_v10_y*.nc')
    #atmo_files = path+'NEMO_ERA5_u10_y2014.nc'

    df = xr.open_mfdataset(atmo_files, parallel=True, chunks={'time':100})

    print(df)

    #get the max for each time step
    max_speed = df.max(dim=["latitude", "longitude"])

    max_speed['v10'].plot()

    plt.title('Max V Velocity')
    plt.tight_layout()

    plt.savefig('ERA5_v10_max.png')

    df.close()

def spatial_max():

    #define a max wind speed, over which we make a spatial plot
    max_value = 60

    figs_path = '/project/6007519/weissgib/plotting/forcing_check_figs/'

    path = '/project/6007519/ANHA4-I/ATMFORCING/ERA5/Ready/'
    u_files = glob.glob(path+'NEMO_ERA5_u10_y*.nc')
    v_files = glob.glob(path+'NEMO_ERA5_v10_y*.nc')

    #u_files = path+'NEMO_ERA5_u10_y2007.nc'
    #v_files = path+'NEMO_ERA5_v10_y2007.nc'

    df_u = xr.open_mfdataset(u_files, parallel=True, chunks={'time':300})
    df_v = xr.open_mfdataset(v_files, parallel=True, chunks={'time':300})

    df = xr.merge([df_u, df_v])
    print(df)

    #calculate the speed
    df['speed'] = np.sqrt(np.square(df['u10'])+np.square(df['v10']))

    #drop all the times where the maximum speed is less than our cutoff
    df_max_u = df.where(df['u10'].max(dim=['latitude', 'longitude']) > max_value, drop=True)
    df_max_v = df.where(df['v10'].max(dim=['latitude', 'longitude']) > max_value, drop=True)

    print(df_max_u)
    print(df_max_v)

    lons = df['longitude']
    lats = df['latitude']
    
    for t in range(df_max_u.dims['time']):
    
        #fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.Robinson()))

        u_max = df_max_u['speed'].isel(time=t)
        print(u_max)

        p1 = u_max.plot(x='longitude', y='latitude', cmap='winter',transform=ccrs.PlateCarree(), subplot_kws=dict(projection=ccrs.Robinson()))
        #p1 = ax.pcolormesh(u_max, lons, lats, transform=ccrs.PlateCarree(), cmap='winter')
        #ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
        #cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')

        #ax.coastlines()
        #plt.title('Speed when U Max > '+str(max_value)+': '+str(t))
        plt.savefig(figs_path+'u_max_spatial_'+str(t)+'.png')
        #plt.show()
        plt.clf()

    df_u.close()
    df_v.close()
    df.close()

if __name__ == "__main__":
    spatial_max()
    #max_timeseries()
