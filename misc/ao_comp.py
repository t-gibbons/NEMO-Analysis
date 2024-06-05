import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#read in the ao index

ao_path = '/mnt/storage6/tahya/obs/monthly_ao_index_b50_current.ascii'

ao = pd.read_csv(ao_path, names=['year', 'month', 'AO index'], header=None, delimiter='   ')

#just want the index between 2002-2019
ao = ao[ao.year >= 2002]
ao = ao[ao.year <= 2020]

#create a datetime column for plotting
ao['day'] = 1

ao['date'] = pd.to_datetime(ao[['year', 'month', 'day']])

print(ao.head)

#take the annual average
mean = ao.groupby(ao.date.dt.year)['AO index'].mean()
print(mean)
