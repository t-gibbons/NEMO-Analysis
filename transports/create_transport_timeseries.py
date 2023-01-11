import glob
import datetime
import scipy.io
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

#get list of all data files

#root_dir = '/project/6000276/weissgib/transports/'
root_dir = '/project/6007519/weissgib/plotting/figs/transports/'
fig_path = '/project/6007519/weissgib/plotting/fwc_figs/'

files = glob.glob(root_dir+'nares_strait_freshwater_transport_*.nc')

#also plot observed values
obs_davis = False 
obs_nares = False
obs_barrow = False

experiment = []
transport = []
date = []

for f in files:
    file_name = f.replace(root_dir, '')
    file_name = file_name.replace('.nc', '')
    exp = file_name[-6:]
    print(exp)
    
    """
    if exp == 'EPM101':
        #exp = 'HYPE,CGRF'
        exp = 'No river water temp'
    if exp == 'EPM102':
        exp = 'HYPE,ERA'
        continue
    if exp == 'EPM014':
        exp = 'Dai and Trenberth,ERA'
        continue
    if exp == 'EPM015':
        exp = 'Dai and Trenberth,CGRF'
        continue
    """
    if exp == 'ETW101':
        continue
        exp = 'River water temp'
   
    d = xr.open_dataset(f)
    #v = (d['vel'].values)
    v = d['__xarray_dataarray_variable__'].values
    datetimeindex = d.indexes['time_counter'].to_datetimeindex()
    times = datetimeindex.values
    l = d.dims['time_counter']

    for i in range(l):
        experiment.append(exp)
    transport.extend(list(v))
    date.extend(list(times))
    d.close()

if obs_davis:
    
    #obs davis strait
    obs_path = '/project/6007519/weissgib/plotting/data_files/observations/DS0413_dOAtrans.mat'

    od = scipy.io.loadmat(obs_path)
    fw_transport = od['FW_dtrans'][:,0]*0.001
    #fw_transport = od['V_dtrans'][:,0]
    dates = od['dates'][0,:]

    timestamps = pd.to_datetime(dates-719529, unit='D')

    transport.extend(list(fw_transport))
    date.extend(list(timestamps))
    for i in range(len(fw_transport)):
        experiment.append('obs')

if obs_nares:
    
    #obs nares strait
    obs_path = '/project/6000276/weissgib/observations/NaresStrait_OBS.mat'

    od = scipy.io.loadmat(obs_path)
    fw_transport = od['fw'][:,0]
    dates = od['date'][:,0]

    timestamps = pd.to_datetime(dates-719529, unit='D')

    transport.extend(list(fw_transport))
    date.extend(list(timestamps))
    for i in range(len(fw_transport)):
        experiment.append('obs')

    """
    #obs nares strait
    obs_path = '/project/6000276/weissgib/observations/Nares-flux2003.dat'
    od = np.genfromtxt(obs_path, comments='%')
    dates = od[:,0]
    fw_transport = od[:,2]*-0.001

    timestamps = pd.to_datetime(dates, unit='D', origin=pd.Timestamp('2003-01-01'))
    print(timestamps)

    transport.extend(list(fw_transport))
    date.extend(list(timestamps))
    for i in range(len(fw_transport)):
        experiment.append('obs')

    obs_path = '/project/6000276/weissgib/observations/Nares-flux2007.dat'
    od = np.genfromtxt(obs_path, comments='%')
    dates = od[:,0]
    fw_transport = od[:,2]*-0.001

    timestamps = pd.to_datetime(dates, unit='D', origin=pd.Timestamp('2003-01-01'))
    print(timestamps)

    transport.extend(list(fw_transport))
    date.extend(list(timestamps))
    for i in range(len(fw_transport)):
        experiment.append('obs')
    """

if obs_barrow:

    #obs barrow strait
    obs_file = '/project/6000276/weissgib/observations/fw_barrow_199808_201008_monthly.txt'
    date_file = '/project/6000276/weissgib/observations/date_barrow_199808_201008_monthly.txt'
    fw_transport = np.loadtxt(obs_file)*-0.000001
    
    t = open(date_file, 'r')
    dates = t.readlines()
    dates = [i.strip() for i in dates]
    t.close()
    timestamps = pd.to_datetime(dates)

    transport.extend(list(fw_transport))
    date.extend(list(timestamps))
    for i in range(len(fw_transport)):
        experiment.append('obs')

#now make a pandas dataframe for easy plotting
all_data = {'experiment': experiment, 'volume_transport': transport, 'date': date}
df = pd.DataFrame(all_data)
print(df)

mean = df.groupby('experiment', as_index=False)['volume_transport'].mean()
print(mean)

#now lets make the time series for each region
rd = df.pivot(index='date', columns='experiment', values='volume_transport')
avg = rd.resample('M').mean() #take the monthly mean

#just output the annual average
annual_avg = rd.resample('Y').mean()
print(annual_avg)

#rd.plot()
#rd["2009-01-03":"2010'12'31"].plot()
avg.plot()
plt.grid(True)
plt.title('Nares Strait')
plt.ylabel('freshwater transport')

#plt.show()
plt.savefig(fig_path+'time_series_freshwater_transport_nares_strait.png')
plt.clf()
