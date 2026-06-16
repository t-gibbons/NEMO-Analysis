import xarray as xr
import matplotlib.pyplot as plt


#here we find the scaling factor
#use the calc clim avg to find a scaling factor for a high and low flow year
#(for Mackenzie, low flow 2019 and high flow 2020)
def calc_scale_factor():

    #read in the arcticGRO data
    obs_files = '/mnt/storage6/tahya/data/ArcticGRO/ArcticGRO_Water_Quality_Data_Mackenzie.csv'

    obs_data = pd.read_csv(obs_files, header=[8])

    #data preprocessing
    obs_data = obs_data[['River', 'Date', 'Discharge', 'Temperature']]
    obs_data = obs_data.iloc[1:]
    obs_data['Date'] = pd.to_datetime(obs_data['Date'], errors='coerce')
    obs_data = obs_data.dropna(subset=['Date'])
    obs_data['Temperature'] = pd.to_numeric(obs_data['Temperature'])
    obs_data['Discharge'] = pd.to_numeric(obs_data['Discharge'])

    print(obs_data)

    obs_data = obs_data.pivot(index='Date', columns='River', values='Temperature')
    print(obs_data)
    clim_mean = obs_data["2003-01-01":"2018-12-31"].mean()
    mean19 = obs_data["2019-01-01":"2019-12-31"].mean()
    mean20 = obs_data["2020-01-01":"2020-12-31"].mean()

    print(clim_mean)
    print(mean19)
    print(mean20)

    scal19 = mean19/clim_mean
    scal20 = mean20/clim_mean

    print(scal19)
    print(scal20)

#scale the 2018 runoff files for a set region
#using the scaling factors that were calculated above
def scale_runoff():

    low_scale=0.75
    high_scale=1.65

    low_year = 2019
    high_year = 2020

    #should just try it based on the region using for camas analysis
    lon_bnds, lat_bnds = (-150, -120), (67, 80)

    runoff_path = '/project/6007519/weissgib/ANHA4-I/RUNOFF/HydroGFD_temp/ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_y2018.nc'

    bf = xr.open_mfdataset(runoff_path, decode_times=False)

    print(bf)

    #now use the bounds to mask the data and only scale in your region
    mask = (bf.nav_lat > lat_bnds[0]) & (bf.nav_lat < lat_bnds[1]) & (bf.nav_lon > lon_bnds[0]) & (bf.nav_lon < lon_bnds[1])

    masked_low = xr.where(mask, low_scale, 1)
    masked_high = xr.where(mask, high_scale, 1)

    low = bf['runoff']*masked_low
    high = bf['runoff']*masked_high

    print(low)

    low.to_netcdf('ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_scaled_y2019.nc')
    high.to_netcdf('ANHA4_ReNat_HydroGFD_HBC_runoff_monthly_scaled_y2020.nc')

    bf.close()


if __name__ == "__main__":
     scale_runoff()

