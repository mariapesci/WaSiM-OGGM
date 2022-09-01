#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 14:06:51 2022

@author: maria

This script contains additional functions required for the coupling
scheme WaSiM-OGGM

"""

import os
import glob
import shutil
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import copy
import os
from netCDF4 import Dataset, date2num
import time
from scipy import stats
import shutil
import tarfile
import gzip
import pickle


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

def clim_adj_anom(rgi_ids, path_in_cru, path_in_user, temp_folder_user, plot=False):
    """
    This function takes the 'climate_historical.nc' file created by OGGM for each glacier
    and using the default climate dataset (i.e. CRU) and uses it as reference for adjusting
    the user's climate dataset through anomalies (i.e. INCA), also using the 
    'climate_historical.nc' files created by OGGM. At the end, new 'climate_historical.nc'
    files are created and saved with the corresponding glacier ID.
    This function is required after creating the glaciers directories and before running any
    climate or calibration functions by OGGM with the user's dataset.
    
    Parameters
    ----------
    rgi_ids : list
        List containing the IDs of the glaciers that are being simulated
    path_in_cru : str
        Path of the mother folder containing the glacier directories from the 
        simulation done with the default climate dataset (i.e. CRU)
    path_in_user : str
        Path of the mother folder containing the glacier directories from the
        simulation done with the user's climate dataset (i.e. INCA)
    temp_folder_user : str
        Path of the folder where the new (adjusted) climate_historical.nc files
        are going to be saved (for all glaciers)
    plot : boolean, optional
        If plots of the adjusted dataset are to be shown for each glacier.
        The default is False.

    Returns
    -------
    The new 'climate_historical.nc' netCDF file for each glacier, with the adjusted
    climate dataset based on anomalies between the default dataset (i.e. CRU) 
    and the user's dataset (i.e. INCA)

    """
    for glac in range(len(rgi_ids)):
        region = rgi_ids[glac].split('.')[0]
        cru_nc = glob.glob(path_in_cru+'/'+region+'/'+region+'.00'+'/'+rgi_ids[glac]+'/climate_historical.nc')[0]
        user_nc = glob.glob(path_in_user+'/'+region+'/'+region+'.00'+'/'+rgi_ids[glac]+'/climate_historical.nc')[0]
                              
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


        x_p_anom = anomalies(meteo_cru2, 'prcp')
        y_p_anom = anomalies(meteo_user2, 'prcp')
        y_p_adj_anom= (((x_p_anom['prcp'] - x_p_anom['mean']) / x_p_anom['std']) * y_p_anom['std'] + y_p_anom['mean']) * x_p_anom['mean']/y_p_anom['mean']


        # Read the netCDF file (original) and adjust the variables 
        ncdata_ref = xr.open_dataset(user_nc)
        ncdata_cru = xr.open_dataset(cru_nc)

        temp_folder_glac = os.path.join(temp_folder_user, rgi_ids[glac])  
        if os.path.exists(temp_folder_glac):
            os.rmdir(temp_folder_glac)
        else:
            os.mkdir(temp_folder_glac)
                
        with Dataset(os.path.join(temp_folder_glac, 'climate_historical.nc') ,'w',format='NETCDF4') as nco: # create the netCDF file
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
        
        if plot is True:
            fig= plt.figure(figsize=(8,12))
            grid = plt.GridSpec(2,1)
            ax1 = plt.subplot(grid[0,0])
            ax2 = plt.subplot(grid[1,0], sharex=ax1)
            ax1.plot(meteo_cru2['temp'].resample('Y').mean(), color='#00008B', label='CRU dataset (OGGM)')
            ax1.plot(meteo_user2['temp'].resample('Y').mean(), color='#FF4500', label='User dataset')
            ax1.plot(y_t_adj_anom.resample('Y').mean(), color='#DC143C', ls='--', label='Adjusted user dataset')
            ax1.set_title('Mean annual temperature values', fontsize=12)
            ax1.set_xlabel('Year', fontsize=12)
            ax1.set_ylabel('Mean annual temperature [°C]', fontsize=12)
            ax1.tick_params(labelsize=12)
 
            ax2.plot(meteo_cru2['prcp'].resample('Y').mean(), color='#00008B', label='CRU dataset (OGGM)')
            ax2.plot(meteo_user2['prcp'].resample('Y').mean(), color='#FF4500', label='User dataset')
            ax2.plot(y_p_adj_anom.resample('Y').mean(), color='#DC143C', ls='--', label='Adjusted user dataset')
            ax2.set_title('Mean annual precipitation values', fontsize=12)
            ax2.set_xlabel('Year', fontsize=12)
            ax2.set_ylabel('Mean annual precipitation [mm]', fontsize=12)
            ax2.tick_params(labelsize=12)
            plt.legend(fontsize=12)
            fig.suptitle('Climate data adjusted through anomalies for glacier '+rgi_ids[glac], fontsize=14)
            plt.show()

def replace_file(path_src, path_dst, rgi_ids, name='climate_historical.nc'):
    """
    This function can be used to update (replace) the climate_historical.nc file
    created 'automatically' by OGGM with the file adjusted after applying the
    anomalies (user's climate dataset adjusted with anomalies according to CRU,
    or any other default dataset).

    Parameters
    ----------
    path_src : str
        path with the location of the climate files that need to be used 
        (user's adjusted dataset)
    path_dst : str
        path with the location of the climate files that need to be replaced
        (created by OGGM, based on CRU or any other default dataset)
    glac_ids: list
        list containing the ID of the glaciers that are simulated

    """
    for glac in range(len(rgi_ids)):
        region = rgi_ids[glac].split('.')[0]
        src = glob.glob(path_src+'/'+rgi_ids[glac]+'/'+name)[0]
        dst = glob.glob(path_dst+'/'+region+'/'+region+'.00'+'/'+rgi_ids[glac]+'/'+name)[0]
        shutil.move(src, dst)
    
    return
