import numpy as np
import numpy.ma as ma
import netCDF4 as nc
from datetime import datetime
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

runids = ['EPM151', 'EPM152', 'EPM014', 'EPM015']

#runids = ['EPM151', 'EPM152']

long_names = {'EPM101': 'Old HYPE, CGRF', 'EPM102': 'Old HYPE, ERA', 'EPM151': 'HYPE, CGRF', 'EPM152': 'HYPE, ERA','EPM014': 'Dai and Trenberth, ERA', 'EPM015': 'Dai and Trenberth, CGRF'}

path = '/project/6007519/weissgib/plotting/data_files/freshwater_content/'
#path = '/project/6000276/weissgib/freshwater_content/'

fig_path = '/project/6007519/weissgib/plotting/figs/timeseries/'

#mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/regions_mask.nc'
mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_temp_regions_mask.nc'
mf = nc.Dataset(mask_file)

experiment = []
region = []
fwc = []
date = []

regions = {'bs_mask': 'Mackenzie Region', 'kara_mask': 'Kara Sea', 'laptev_mask': 'Laptev Sea', 'bs_east_mask': 'Eastern Bering Strait'}

#regions = {'caa_mask': 'Canadian Arctic Archipelago', 'ca_mask': 'Central Arctic', 'cs_mask': 'Canadian Shelf', 'cb_mask': 'Canadian Basin', 'eb_mask': 'Eurasian Basin', 'ss_mask': 'Siberian Shelf', 'ds_mask': 'Davis Strait', 'hb_mask': 'Hudson Bay', 'bs_mask': 'Bering Strait', 'ls_mask': 'Labrador Sea', 'ns_mask': 'Nares Strait', 'fs_mask': 'Fram Strait', 'lc_mask': 'Labrador Current'}

for exp in runids:
    fwc_file = path+exp+'_fwc_34.8_isohaline_0m.nc'
    
    df = xr.open_mfdataset(fwc_file)
    
    fwc_data = df['vosaline'].values
    times = df['time_counter'].values
    t = df.dims['time_counter']
    
    df.close()

    times = [datetime.strptime(str(x),'%Y-%m-%d %H:%M:%S' ) for x in list(times)]
    
    for r in regions.keys():
        rmask = mf[r][0,0]

        #get the mask into the correct shape
        rmask = np.broadcast_to(rmask,(t,)+rmask.shape)

        masked_fwc = ma.masked_where(rmask==1,fwc_data)
        regional_fwc_ts = masked_fwc.mean(axis=(1,2))
        
        for i in range(t):
            experiment.append(long_names[exp])
            region.append(r)
        
        fwc.extend(list(regional_fwc_ts))
        date.extend(times)
        

#make a pandas database for plotting
all_data = {'experiment': experiment, 'region':region, 'fwc':fwc, 'date':date}
df = pd.DataFrame(all_data)

#now we can do comparison timeseries plots
for r in regions:
    rd = df.loc[df['region'] == r]
    rd = rd.pivot(index='date', columns='experiment', values='fwc')
    rd.plot()
    plt.grid(True)
    plt.title(regions[r]+' Freshwater Content')
    plt.ylabel('freshwater content (m)')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
    plt.legend()
    plt.tight_layout()
    #plt.show()
    plt.savefig(fig_path+r+'_fwc_34.8_isohaline_timeseries_new.png')
    plt.clf()
