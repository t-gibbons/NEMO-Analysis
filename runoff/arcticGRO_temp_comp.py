"""
author: Tahya Weiss-Gibbons, weissgib@ualberta.ca

Plots the average temperature for each month for the 6 major rivers
Observational data from the ArcticGRO database
Uses raw A-HYPE model ouput, taken at the same locations as the sampling from ArcticGRO
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist


#first read in all of the observational files
obs_path = '/mnt/storage6/tahya/obs/ArcticGRO/'

obs_files = glob.glob(os.path.join(obs_path, "ArcticGRO*.csv"))

obs_data = pd.concat((pd.read_csv(f, header=[8]) for f in obs_files), ignore_index=True)

#data preprocessing
obs_data = obs_data[['River', 'Date', 'Discharge', 'Temperature']]
obs_data = obs_data.iloc[1:]
obs_data['Date'] = pd.to_datetime(obs_data['Date'], errors='coerce')
obs_data = obs_data.dropna(subset=['Date'])
obs_data['Temperature'] = pd.to_numeric(obs_data['Temperature'])

print(obs_data)

#now lets get the raw ahype data
hype_path = '/mnt/storage4/tahya/runoff/new_runoff_remapping/HYDRO-GFD/NEMO_HYPE_Results_v2_Arctic_Calib_temp_1981_2019_HydroGFD.txt'

hype_data = pd.read_csv(hype_path, sep='\t')

#data preprocessing
hype_data = hype_data.T
hype_data.columns = hype_data.iloc[0]
hype_data['SUBID'] = hype_data.index

hype_data = hype_data.drop(hype_data.index[:1])
hype_data = hype_data.apply(pd.to_numeric)

#remove data from before the arcticGRO period
old_dates = hype_data.filter(regex='19.*', axis=1)
hype_data = hype_data.drop(old_dates.columns.values, axis=1)
print(hype_data)

#river locations are taken from the ArcticGRO metadata
river_data = {'River Name': ['Kolyma', 'Lena', 'Mackenzie', "Ob'", 'Yensiey', 'Yukon'], 'Lat': [68.75, 66.77, 67.45, 66.63, 69.38, 61.93], 'Lon': [161.30, 123.37, 133.74, 66.60, 86.15, 162.88]}
river_loc = pd.DataFrame(data=river_data)

#now lets find the closest lat lon location from a-hype
distance = cdist(river_loc[["Lat", "Lon"]], hype_data[["LAT", "LON"]], metric="euclidean")

river_loc['Closest A-Hype'] = hype_data['SUBID'].to_numpy()[distance.argmin(axis=1)]

print(river_loc)

#average everything monthly
obs_monthly = obs_data.groupby([obs_data['Date'].dt.month, "River"])['Temperature'].mean().reset_index()

print(obs_monthly)

#plot each rivers seasonal cycle against the a-hype average
for index, row in river_loc.iterrows():
    river = row['River Name']
    print(river)
    hype_river = hype_data.loc[hype_data['SUBID'] == row['Closest A-Hype']]

    hype_river = hype_river.T
    hype_river = hype_river.iloc[3:]
    hype_river = hype_river.iloc[:-1]
    hype_river.index = pd.to_datetime(hype_river.index)

    hype_monthly = hype_river.groupby([hype_river.index.month]).mean().reset_index()
    hype_monthly = hype_monthly.rename({'SUBID': 'Date', hype_monthly.columns[1]: 'Temperature'}, axis=1)
    hype_monthly[hype_monthly < 0] = 0

    obs_river = obs_monthly.loc[obs_monthly['River'] == river]

    print(hype_monthly)
    print(obs_river)

    plt.plot(hype_monthly['Date'], hype_monthly['Temperature'], label='A-Hype')
    plt.plot(obs_river['Date'], obs_river['Temperature'], label='ArcticGRO')
    plt.title(river+" average seasonal temperature")
    plt.legend()

    plt.savefig(river+'_temp_seasonal.png')
    plt.clf()
