import glob
import datetime
import scipy.io
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

#get list of all data files

root_dir = '/project/6007519/weissgib/plotting/figs/transports/'
fig_path = '/project/6007519/weissgib/plotting/transports/'

files = glob.glob(root_dir+'davis_strait_freshwater_transport_*.nc')

#also plot observed values
obs_davis = True
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
    
    
    if exp == 'EPM101':
        continue
        #exp = 'HYPE,CGRF'
        exp = 'No river water temp'
        continue
    if exp == 'EPM102':
        exp = 'HYPE,ERA'
        continue
    if exp == 'EPM151':
        #exp = 'No river water temp'
        exp = 'A-HYPE'
    if exp == 'EPM152':
        continue
        #exp = 'HYPE, ERA'
    if exp == 'EPM014':
        continue
        #exp = 'Dai and Trenberth,ERA'
    if exp == 'EPM015':
        #continue
        exp = 'Dai and Trenberth'
    if exp == 'ETW161':
        exp = 'River water temp'
        continue
    
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
        experiment.append('Observations')

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

#df['diff'] = df['River water temp']-df['No river water temp']

mean = df.groupby('experiment', as_index=False)['volume_transport'].mean()
print(mean)

#now lets make the time series for each region
rd = df.pivot(index='date', columns='experiment', values='volume_transport')
avg = rd.resample('M').mean() #take the monthly mean
#avg = avg.loc["2008-01-01":"2018-12-31"]
print(avg)
avg['diff'] = ((abs(avg['A-HYPE'])-abs(avg['Dai and Trenberth']))/abs(avg['Dai and Trenberth']))*100
mean = avg['diff'].mean()
print(mean)
#mean15 = avg['Dai and Trenberth'].mean()
#mean151 = avg['A-HYPE'].mean()
#meanobs = avg['Observations'].mean()
#print(mean15)
#print(mean151)
#print(meanobs)

#just output the annual average
annual_avg = rd.resample('Y').mean()
print(annual_avg)

#rd.plot()
#rd["2009-01-03":"2010'12'31"].plot()
avg['diff'].plot()
plt.axhline(y=mean, color='C0', linestyle='--')
#avg.plot()
#plt.axhline(y=mean151, color='C0', linestyle='--')
#plt.axhline(y=mean15, color='C1', linestyle='--')
#plt.axhline(y=meanobs, color='C2', linestyle='--')
ax = plt.gca()
ax.set_ylim([-30, 30])
plt.grid(True)
plt.title('Davis Strait')
plt.ylabel('percentage change freshwater transport')

#plt.show()
plt.savefig(fig_path+'time_series_freshwater_transport_davis_strait_CGRF_change.png')
plt.clf()
