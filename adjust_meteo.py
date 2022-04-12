#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 

@author: María Pesci
"""

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
from scipy import stats

def extract_time_series(nc_file_path, var_name, multiple=False, time_bnds=False):
    df =  pd.DataFrame()
    if multiple:
        if time_bnds:
            ncdata = xr.open_mfdataset(nc_file_path+'/'+var_name+'*.nc', drop_variables='time_bnds')
        else:
            ncdata = xr.open_mfdataset(nc_file_path+'/'+var_name+'*.nc')
        index = ncdata.coords['time'].to_dataframe()['time']
        var = ncdata.variables[var_name]
    else:
        ncdata = xr.open_dataset(nc_file_path)
        index = ncdata.coords['time'].to_dataframe()['time']
        var = ncdata.variables[var_name]

    df = pd.DataFrame(index=index)
    df[str(var_name)] = var[:].values

    return df

# ===== INPUTS & DEFINITIONS =====
# Path containing the "default" climate" data from OGGM (i.e. CRU)
path_in_cru = os.path.abspath(r'E:\DIRT_X\OGGM\meteo_reg\climate_historical_cru.nc')
# Path containing the 
path_in_user = os.path.abspath(r'E:\DIRT_X\OGGM\meteo_reg\climate_historical_inca.nc') ##### CHANGE THE PATH #####
#path_in_user = os.path.abspath(r'E:\DIRT_X\OGGM\meteo_reg\climate_historical_st.nc') # for stations
#path_in_user = os.path.abspath(r'E:\DIRT_X\OGGM\Gepatschferner_nc_files\INCA\INCA_wasim_normal\climate_historical.nc') # for INCA
                               
# For Vernagtferner!!!
# path_in_cru = os.path.abspath(r'E:\DIRT_X\OGGM\Vernagtferner\CRU\climate_historical.nc') ##### CHANGE THE PATH #####
# path_in_user = os.path.abspath(r'E:\DIRT_X\OGGM\Vernagtferner\INCA\climate_historical.nc') ##### CHANGE THE PATH #####

# FOR CRU   
prcp_cru = extract_time_series(path_in_cru, 'prcp')
temp_cru = extract_time_series(path_in_cru, 'temp')
meteo_cru = pd.concat([prcp_cru, temp_cru], axis=1)

# FOR INCA   
prcp_user = extract_time_series(path_in_user, 'prcp')
temp_user = extract_time_series(path_in_user, 'temp')
meteo_user = pd.concat([prcp_user, temp_user], axis=1)

t1 = max(prcp_cru.index[0], prcp_user.index[0])
t2 = min(prcp_cru.index[-1], prcp_user.index[-1])

meteo_cru2 =  meteo_cru[(meteo_cru.index >= t1) & (meteo_cru.index <= t2)]
meteo_user2 =  meteo_user[(meteo_user.index >= t1) & (meteo_user.index <= t2)]

# TEMPERATURE
x_t = meteo_cru2['temp'].to_numpy()
y_t = meteo_user2['temp'].to_numpy()
slope_t, intercept_t, r_value_t, p_value_t, std_err_t = stats.linregress(x_t, y_t)

# Plot outputs
line_t = slope_t*x_t+intercept_t
fig, ax = plt.subplots(figsize=(12,12))
plt.scatter(x_t, y_t,  color='#808080')
plt.plot(x_t, line_t, 'r', label='y={:.2f}x+{:.2f}'.format(slope_t,intercept_t))
plt.title('Linear regression for temperature', fontsize='xx-large')
plt.xlabel('Mean monthly temperature CRU - Observed [°C]', fontsize='x-large')
plt.ylabel('Mean monthly temperature User - simulated [°C]', fontsize='x-large')
plt.legend(fontsize=16)
plt.xlim(-17,10)
plt.ylim(-17,10)
plt.show()

# Adjustment of the user dataset (x_t_adj) to the default (CRU) dataset:
x_t_adj = (y_t - intercept_t) / slope_t
meteo_user2['temp_adj'] = x_t_adj

fig, ax = plt.subplots(figsize=(12,12))
ax.plot(meteo_cru2['temp'].resample('Y').mean(), color='#00008B', label='CRU dataset (OGGM)')
ax.plot(meteo_user2['temp'].resample('Y').mean(), color='#FF4500', label='User dataset')
ax.plot(meteo_user2['temp_adj'].resample('Y').mean(), color='#DC143C', ls='--', label='Adjusted user dataset')
plt.title('Mean annual temperature values', fontsize='xx-large')
plt.xlabel('Year', fontsize='x-large')
plt.ylabel('Mean annual temperature [°C]', fontsize='x-large')
plt.legend(fontsize=14)
plt.show()

# PRECIPITATION
t1_p = t1 #'2016-01-01 00:00:00'
t2_p = t2
x_p = meteo_cru2['prcp'][t1_p:t2_p].to_numpy()
y_p = meteo_user2['prcp'][t1_p:t2_p].to_numpy()
slope_p, intercept_p, r_value_p, p_value_p, std_err_p = stats.linregress(x_p, y_p)

# Simple multiplicative correction
alfa = y_p.mean() / x_p.mean()

# Adjustment of the user dataset (x_p_adj) to the default (CRU) dataset:
x_p_adj = y_p / alfa

# Plot outputs
line_p = slope_p*x_p+intercept_p
fig, ax = plt.subplots(figsize=(12,12))
plt.scatter(x_p, y_p,  color='black')
plt.plot(x_p, line_p, 'r', label='y={:.2f}x+{:.2f}'.format(slope_p,intercept_p))
plt.title('Linear regression for precipitation', fontsize='x-large')
plt.xlabel('Mean monthly precipitation CRU - Observed [mm]', fontsize='x-large')
plt.ylabel('Mean monthly precipitation User - simulated [mm]', fontsize='x-large')
plt.xlim(0,450)
plt.ylim(0,450)
#plt.legend(fontsize=12)
plt.show()

meteo_user2['prcp_adj'] = x_p_adj

fig, ax = plt.subplots(figsize=(12,12))
ax.plot(meteo_cru2['prcp'].resample('Y').mean(), color='#00008B', label='CRU dataset (OGGM)')
ax.plot(meteo_user2['prcp'].resample('Y').mean(), color='#90EE90', label='User dataset')
ax.plot(meteo_user2['prcp_adj'].resample('Y').mean(), color='#15B01A', ls='--', label='Adjusted user dataset')
plt.title('Mean annual precipitation values', fontsize='xx-large')
plt.xlabel('Year', fontsize='x-large')
plt.ylabel('Mean annual precipitation [mm]', fontsize='x-large')
plt.legend(fontsize=14)
plt.show()