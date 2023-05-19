"""
fwc_region_timeseries.py
Author: Tahya Weiss-Gibbons (weissgib@ualberta.ca)

plot the freshwater content in a region as a timeseries
can either plot a number of model runs against each other on same plot
or the difference between two model runs
"""

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4 as nc
import numpy.ma as ma
from datetime import datetime
import matplotlib.pyplot as plt

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

#define variables used for all functions

long_names = {'EPM101': 'Old HYPE, CGRF', 'EPM102': 'Old HYPE, ERA', 'EPM151': 'HYPE, CGRF', 'EPM152': 'HYPE, ERA','EPM014': 'Dai and Trenberth, ERA', 'EPM015': 'Dai and Trenberth, CGRF'}

path = '/project/6007519/weissgib/plotting/data_files/freshwater_content/'
fig_path = '/project/6007519/weissgib/plotting/fwc_figs/'
#mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/runoff_temp_regions_mask.nc'
mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/regions_mask.nc'

#regions = {'bs_mask': 'Mackenzie Region', 'kara_mask': 'Kara Sea', 'laptev_mask': 'Laptev Sea', 'bs_east_mask': 'Eastern Bering Strait'}
regions = {'caa_mask': 'Canadian Arctic Archipelago', 'ca_mask': 'Central Arctic', 'cs_mask': 'Canadian Shelf', 'cb_mask': 'Canadian Basin', 'eb_mask': 'Eurasian Basin', 'ss_mask': 'Siberian Shelf', 'ds_mask': 'Davis Strait', 'hb_mask': 'Hudson Bay', 'bs_mask': 'Bering Strait', 'ls_mask': 'Labrador Sea', 'ns_mask': 'Nares Strait', 'fs_mask': 'Fram Strait', 'lc_mask': 'Labrador Current'}

def read_mask(runids, long_name):

    mf = nc.Dataset(mask_file)

    experiment = []
    region = []
    fwc = []
    date = []

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
            print(regional_fwc_ts.shape)
        
            for i in range(t):

                if long_name: 
                    experiment.append(long_names[exp])
                else: 
                    experiment.append(exp)
                region.append(r)
        
            fwc.extend(list(regional_fwc_ts))
            date.extend(times)
        

    #make a pandas database for plotting
    all_data = {'experiment': experiment, 'region':region, 'fwc':fwc, 'date':date}
    df = pd.DataFrame(all_data)

    return df

def fwc_region_comp(runids, long_name=False):

    df = read_mask(runids, long_name)

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

def fwc_region_diff(runids, long_name=False):

    df = read_mask(runids, long_name)

    #plot the difference between two runs
    for r in regions:
        rd = df.loc[df['region'] == r]
        rd = rd.pivot(index='date', columns='experiment', values='fwc')
        rd['diff'] = rd[runids[0]] - rd[runids[1]]
        rd['diff'].plot()
        plt.grid(True)
        plt.title(regions[r]+' Freshwater Content Difference')
        plt.ylabel('freshwater content (m)')
        #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
        plt.legend()
        plt.tight_layout()
        #plt.show()
        plt.savefig(fig_path+r+'_fwc_34.8_isohaline_timeseries_diff_cgrf.png')
        plt.clf()


if __name__ == "__main__":

    fwc_region_diff(runids = ['EPM151', 'EPM015'])
