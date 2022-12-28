"""
Using the calculated freshwater from a model run, compares with the observations for the beaufort gyre
Uses the BG definition from Proshutinsky 
"""
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import numpy.ma as ma
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

###this is for running script backend###
import matplotlib
matplotlib.use('Agg')
###----------------------------------###

path = '/project/6007519/weissgib/plotting/data_files/freshwater_content/'
fig_path = '/project/6007519/weissgib/plotting/fwc_figs/'
mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/regions_mask.nc'
grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'

runids = ['EPM151', 'EPM152', 'EPM014', 'EPM015']

mesh = nc.Dataset(grid_file)

e1t = np.array(mesh.variables['e1t'])[0,:,:]
e2t = np.array(mesh.variables['e2t'])[0,:,:]

mesh.close()

mask = nc.Dataset(mask_file)

cb_mask = mask['cb_mask'][0,0]

mask.close()

dates = []
fwc_all = []
experiment = []


for r in runids:
    model_file = path+r+'_fwc_34.8_isohaline_0m.nc'
    
    mf = xr.open_mfdataset(model_file)
    
    x = mf.dims['x_grid_T']
    y = mf.dims['y_grid_T']
    t = mf.dims['time_counter']
   
    #get the mask into the correct shape
    mask = np.broadcast_to(cb_mask,(t,)+cb_mask.shape)

    times = mf['time_counter'].values
    times = [(datetime.strptime(str(x),'%Y-%m-%d %H:%M:%S' )).date() for x in list(times)]
    dates.extend(times)
    
    for i in range(t):
        experiment.append(r)
    
    fwc = mf['vosaline'].values
    fwc_area = np.zeros((t,y,x))
    
    for j in range(y):
        for i in range(x):
            fwc_area[:,j,i] = fwc[:,j,i]*e1t[j,i]*e2t[j,i]*10**(-12)
    
    #now just average over the canadian basin area
    masked_fwc = ma.masked_where(mask==1, fwc_area)

    fwc_ts = masked_fwc.sum(axis=(1,2)) #should this be the mean or the sum over the area??
    
    fwc_all.extend(list(fwc_ts))
    
    mf.close()


#and read in the observational data
obs_file = '/project/6007519/weissgib/plotting/data_files/observations/beaufort_gyre_obs_proshutinsky_2019_ssh_derived.csv'
of = pd.read_csv(obs_file, delim_whitespace=True)
years = of['Year'].tolist()
date = [(datetime(int(x),1,1) + timedelta(days=(x%1)*365)).date() for x in years]

dates.extend(date)

for i in range(len(date)):
    experiment.append('obs')

fwc_all.extend(of['FWC1'].tolist())

all_data = {'experiment': experiment, 'date': dates, 'fwc': fwc_all}
df = pd.DataFrame(all_data)

print(df)

df['date'] = pd.to_datetime(df['date'])

rd = df.pivot(index='date', columns='experiment', values='fwc')

avg = rd.resample('M').mean()

plt.plot(rd['obs'], label='observations')

for r in runids:
    plt.plot(rd[r], label=r)

plt.legend()

#avg.plot(x_compat=True)
plt.grid(True)
plt.title('Beaufort Gyre Freshwater Content')
plt.tight_layout()
plt.savefig(fig_path+'cb_fwc_timeseries_new_obs.png')


