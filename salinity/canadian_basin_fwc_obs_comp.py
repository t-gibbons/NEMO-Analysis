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


#set your paths for inputs and figure output here
path = '/project/6007519/weissgib/plotting/data_files/freshwater_content/'
fig_path = '/project/6007519/weissgib/plotting/fwc_figs/'
mask_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/regions_mask.nc'
grid_file = '/project/6007519/weissgib/plotting/data_files/anha4_files/ANHA4_mesh_mask.nc'

#runs which you have the freshwater content calculated that you want to compare
runids = ['EPM151', 'EPM152', 'EPM014', 'EPM015']

exp_names = {'EPM151': 'HYPE, CGRF', 'EPM152': 'HYPE, ERA', 'EPM014': 'Dai, ERA', 'EPM015': 'Dai, CGRF'} #long names

#need information about the model grid
mesh = nc.Dataset(grid_file)

e1t = np.array(mesh.variables['e1t'])[0,:,:]
e2t = np.array(mesh.variables['e2t'])[0,:,:]

mesh.close()

#and need the bg mask
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
        experiment.append(exp_names[r])
    
    fwc = mf['vosaline'].values
    fwc_area = np.zeros((t,y,x))
    
    for j in range(y):
        for i in range(x):
            fwc_area[:,j,i] = fwc[:,j,i]*e1t[j,i]*e2t[j,i]*10**(-12)
    
    #now just average over the canadian basin area
    masked_fwc = ma.masked_where(mask==1, fwc_area)

    fwc_ts = masked_fwc.sum(axis=(1,2))
    
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

#add all data to a pandas dataframe for plotting
all_data = {'experiment': experiment, 'date': dates, 'fwc': fwc_all}
df = pd.DataFrame(all_data)

df['date'] = pd.to_datetime(df['date'])

rd = df.pivot(index='date', columns='experiment', values='fwc')

#need to make sure the observations are on the same date spacing so plots correctly
avg = rd.resample('M').mean()
rd['obs'] = rd['obs'].interpolate(method='linear', limit=8)

plt.plot(rd['obs'], label='observations')

for r in runids:
    plt.plot(rd[exp_names[r]], label=exp_names[r])

plt.legend()

#avg.plot(x_compat=True)
plt.grid(True)
plt.title('Beaufort Gyre Freshwater Content')
plt.tight_layout()
plt.savefig(fig_path+'cb_fwc_timeseries_new_obs.png')


