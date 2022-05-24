#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 

@author: María Pesci

This script is the same as meteo_reg.py but adapted to work on a list of glaciers
(for example, for Gepatschalm with 16 glaciers).

"""

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
from netCDF4 import Dataset, date2num
import time
from scipy import stats
import glob

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

def anomalies(x,variable):
    '''
    Takes a time series of monthly data (temperature or precipitation) from 
    a reference dataset (i.e. CRU) or a user's dataset and calculates the
    long-term mean and corresponding empirical standard deviation.

    Parameters
    ----------
    x : time series of meteorological data (i.e. CRU)
    variable : str, precipitation or temperature

    Returns
    -------
    mean and standard deviation values (for anomalies) of a given dataset
    '''
    
    x_mean = (x[variable].groupby('{:%m-%d}'.format).mean()).to_dict()
    x_std = (x[variable].groupby('{:%m-%d}'.format).std()).to_dict()
    x_o = pd.DataFrame(copy.deepcopy(x))
    x_o['month_day'] = x_o.index.strftime('%m-%d')
    x_o['mean'] = 0.
    x_o['std'] = 0.
    for i in range(len(x_o['month_day'])):
        if x_o['month_day'][i] in x_mean:
           x_o['mean'][i] = x_mean[x_o['month_day'][i]]
           x_o['std'][i] = x_std[x_o['month_day'][i]]
    
    return x_o

# definitions
path_in_cru = os.path.abspath(r'E:\DIRT_X\OGGM\Gepatschalm_nc_files\CRU') ##### CHANGE THE PATH #####
path_in_user = os.path.abspath(r'E:\DIRT_X\OGGM\Gepatschalm_nc_files\INCA\INCA_wasim_normal') ##### CHANGE THE PATH #####

rgi_ids = ['RGI60-11.00746'] #, 'RGI60-11.00770', 'RGI60-11.00783', 'RGI60-11.00768', 
          # 'RGI60-11.00732', 'RGI60-11.00725', 'RGI60-11.00698', 'RGI60-11.00708',
          # 'RGI60-11.00717', 'RGI60-11.00714', 'RGI60-11.00711', 'RGI60-11.00684', 
          # 'RGI60-11.00712', 'RGI60-11.00756', 'RGI60-11.00716', 'RGI60-11.00709']

for glac in range(len(rgi_ids)):
    cru_nc = glob.glob(path_in_cru+'\\'+rgi_ids[glac]+'\climate_historical.nc')[0]
    user_nc = glob.glob(path_in_user+'\\'+rgi_ids[glac]+'\climate_historical.nc')[0]
                          
    # FOR CRU   
    prcp_cru = extract_time_series(cru_nc, 'prcp')
    temp_cru = extract_time_series(cru_nc, 'temp')
    meteo_cru = pd.concat([prcp_cru, temp_cru], axis=1)
    
    # FOR USER  
    prcp_user = extract_time_series(user_nc, 'prcp')
    temp_user = extract_time_series(user_nc, 'temp')
    meteo_user = pd.concat([prcp_user, temp_user], axis=1)

    t1 = max(prcp_cru.index[0], prcp_user.index[0])
    t2 = min(prcp_cru.index[-1], prcp_user.index[-1])

    meteo_cru2 =  meteo_cru[(meteo_cru.index >= t1) & (meteo_cru.index <= t2)]
    meteo_user2 =  meteo_user[(meteo_user.index >= t1) & (meteo_user.index <= t2)]


    # Using anomalies (Förster et al., 2018)
    y_obs = meteo_cru2['temp']
    x_t_anom = anomalies(meteo_cru2, 'temp')
    y_t_anom = anomalies(meteo_user2, 'temp')
    y_t_adj_anom= (((x_t_anom['temp'] - x_t_anom['mean']) / x_t_anom['std']) * y_t_anom['std'] + y_t_anom['mean']) * x_t_anom['mean']/y_t_anom['mean']

    # Plot outputs
    fig, ax = plt.subplots(figsize=(12,12))
    ax.plot(meteo_cru2['temp'].resample('Y').mean(), color='#00008B', label='CRU dataset (OGGM)')
    ax.plot(meteo_user2['temp'].resample('Y').mean(), color='#FF4500', label='User dataset')
    ax.plot(y_t_adj_anom.resample('Y').mean(), color='#DC143C', ls='--', label='Adjusted user dataset')
    plt.title('Mean annual temperature values', fontsize=20)
    plt.xlabel('Year', fontsize=20)
    plt.xticks(fontsize=20)
    plt.ylabel('Mean annual temperature [°C]', fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=20)
    #plt.savefig(r'E:\Workshops_presentations\Montpellier_2022\temp_adj.png', dpi=1000)
    plt.show()
    
    x_p_anom = anomalies(meteo_cru2, 'prcp')
    y_p_anom = anomalies(meteo_user2, 'prcp')
    y_p_adj_anom= (((x_p_anom['prcp'] - x_p_anom['mean']) / x_p_anom['std']) * y_p_anom['std'] + y_p_anom['mean']) * x_p_anom['mean']/y_p_anom['mean']
    
    # Plot outputs
    fig, ax = plt.subplots(figsize=(12,12))
    ax.plot(meteo_cru2['prcp'].resample('Y').mean(), color='#00008B', label='CRU dataset (OGGM)')
    ax.plot(meteo_user2['prcp'].resample('Y').mean(), color='#FF4500', label='User dataset')
    ax.plot(y_p_adj_anom.resample('Y').mean(), color='#DC143C', ls='--', label='Adjusted user dataset')
    plt.title('Mean annual precipitation values', fontsize=20)
    plt.xlabel('Year', fontsize=20)
    plt.xticks(fontsize=20)
    plt.ylabel('Mean annual precipitation [mm]', fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=20)
    #plt.savefig(r'E:\Workshops_presentations\Montpellier_2022\prec_adj.png', dpi=1000)
    plt.show()

    # Read the netCDF file (original) and adjust the variables 
    ncdata_ref = xr.open_dataset(path_in_user+'\\'+rgi_ids[glac]+'\climate_historical.nc')
    ncdata_cru = xr.open_dataset(path_in_cru+'\\'+rgi_ids[glac]+'\climate_historical.nc')

    with Dataset(os.path.join(r'E:\DIRT_X\OGGM\meteo_reg\Gepatschalm\climate_historical_anom'+rgi_ids[glac]+'.nc') ,'w',format='NETCDF4') as nco: # create the netCDF file
        nco.ref_hgt = ncdata_cru.attrs['ref_hgt']
        nco.ref_pix_lon = ncdata_cru.attrs['ref_pix_lon']
        nco.ref_pix_lat = ncdata_cru.attrs['ref_pix_lat']
        nco.ref_pix_dis = ncdata_cru.attrs['ref_pix_dis']
        nco.climate_source = ncdata_ref.attrs['climate_source']
        nco.hydro_yr_0 = ncdata_cru.attrs['hydro_yr_0']
        nco.hydro_yr_1 = ncdata_cru.attrs['hydro_yr_1']
        nco.author = ncdata_ref.attrs['author']
        nco.author_info = ncdata_ref.attrs['author_info']
        
        # Add dimension
        dim_time = nco.createDimension('time', None) # in this case the time dimension is "unlimited"
        # Add variables
        
        nc_time = nco.createVariable('time', 'i4', ('time'))
        nc_time.units = 'days since 1801-01-01 00:00:00'
        nc_time.calendar = 'standard'
        time = meteo_cru2.index
        numdate = date2num([t for t in time], nc_time.units, calendar=nc_time.calendar)
        nc_time[:] = numdate
        
        # create variable for precipitation data
        var_prcp = nco.createVariable('prcp', 'f4',  ('time'))
        var_prcp.units = 'kg m-2'
        var_prcp.long_name = 'total monthly precipitation amount'
        var_prcp[:] = y_p_adj_anom #x_p_adj #
        
        # create variable for temperature data
        var_temp = nco.createVariable('temp', 'f4',  ('time'))
        var_temp.units = 'degC'
        var_temp.long_name = '2m temperature at height ref_hgt'
        var_temp[:] = y_t_adj_anom # meteo_user2['temp_adj'] #

