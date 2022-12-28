import datetime
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt

def bio_var_region(variable):

    runids = ['EPM101', 'EPM102', 'EPM014', 'EPM015']
    path = "/project/6007519/weissgib/plotting/data_files/surf_bio/"
    fig_path = '/project/6007519/weissgib/plotting/figs/bio_var/'

    regions = {'caa_mask': 'CAA'}
    #regions = {'caa_mask': 'Canadian Arctic Archipelago', 'ca_mask': 'Central Arctic', 'cs_mask': 'Canadian Shelf', 'cb_mask': 'Canadian Basin', 'eb_mask': 'Eurasian Basin', 'ss_mask': 'Siberian Shelf', 'ds_mask': 'Davis Strait', 'hb_mask': 'Hudson Bay', 'bs_mask': 'Bering Strait', 'ls_mask': 'Labrador Sea', 'ns_mask': 'Nares Strait', 'fs_mask': 'Fram Strait', 'lc_mask': 'Labrador Current'}

    long_name = {'vooxy': 'dissolved oxygen', 'vodic': 'dissolved carbon dioxide', 'voalk': 'alkalinity', 'pco2': 'atmos-ocean carbon dioxide flux'}

    experiment = []
    date = []
    var = []
    region = []

    for r in regions.keys():
        for runid in runids:
            d = xr.open_mfdataset(path+runid+'_'+variable+'_'+r+'.nc')

            v = d[variable].values
            #v = d['pco2_surf'].values

            datetimeindex = d.indexes['time_counter'].to_datetimeindex()
            times = datetimeindex.values
            l = d.dims['time_counter']

            if runid == 'EPM101':
                runid = 'HYPE,CGRF'
            if runid == 'EPM102':
                runid = 'HYPE,ERA'
                continue
            if runid == 'EPM014':
                runid = 'Dai and Trenberth,ERA'
                continue
            if runid == 'EPM015':
                runid = 'Dai and Trenberth,CGRF'

            for i in range(l):
                region.append(r)
                experiment.append(runid)
            var.extend(list(v))
            date.extend(list(times))

            d.close()

    #now make a pandas dataframe for easy plotting
    all_data = {'experiment': experiment, 'region': region, variable: var, 'date': date}
    df = pd.DataFrame(all_data)
    print(df)

    #now lets make the time series for each region
    for r in regions:
        rd = df.loc[df['region'] == r]
        rd = rd.pivot(index='date', columns='experiment', values=variable)
        rd.plot()
        plt.grid(True)
        plt.title(regions[r]+' '+long_name[variable])
        plt.ylabel('mol/m^3')
        #plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=2, fancybox=True, shadow=True)
        plt.legend()
        plt.tight_layout()
        #plt.show()
        plt.savefig(fig_path+variable+'_'+r+'_cgrf_timeseries.png')
        plt.clf()

if __name__ == "__main__":
    bio_var_region('vooxy')
    bio_var_region('vodic')
    bio_var_region('voalk')
    #bio_var_region('pco2')

