import glob
import datetime
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

root_dir = '/mnt/storage6/tahya/model_files/'

section = 'davis_strait'

runids = ['EPM015', 'EPM151']

#want the freshwater and volume flux

experiment = []
freshwater = []
volume = []
date = []

for r in runids:
    ff = root_dir+section+'_freshwater_transport_'+r+'.nc'
    vf = root_dir+section+'_volume_transport_'+r+'.nc'

    dv = xr.open_dataset(vf)
    v = dv['vel'].values
    datetimeindex = dv.indexes['time_counter'].to_datetimeindex()
    times = datetimeindex.values
    l = dv.dims['time_counter']

    dv.close()

    df = xr.open_dataset(ff)
    f = df['__xarray_dataarray_variable__'].values

    df.close()
    
    for i in range(l):
        experiment.append(r)

    volume.extend(list(v))
    freshwater.extend(list(f))
    date.extend(list(times))

all_data = {'experiment':experiment, 'volume_transport':volume, 'freshwater_transport':freshwater,'date':date}

dataframe = pd.DataFrame(all_data)
print(dataframe)

df = dataframe.pivot(index='date', columns='experiment')
avg = df.resample('M').mean()
print(avg)

fig, ax = plt.subplots()

avg.plot(y='volume_transport', ax=ax, color=['C0','C1'])
avg.plot(y='freshwater_transport', ax=ax, secondary_y=True, linestyle='dashed', color=['C0','C1'])

ax.set_ylabel('Volume Transport (Sv)')
ax.right_ax.set_ylabel('Freshwater Transport (Sv)')

ax.legend()
plt.title('Davis Strait Transport')
plt.tight_layout()

#plt.show()
plt.savefig(section+'_vol_fresh_timeseries.png')

  
