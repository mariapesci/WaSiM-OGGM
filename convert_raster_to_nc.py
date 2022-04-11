#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 

@author: María Pesci

This script reads the monthly output from WaSiM for temperature and precipitation, which are provided in ASCII format
and writes one netCDF file containing both variables, to be read by OGGM.
Moreover, the coordinates and altitude of each grid cell is read and assigned to the files.
As an addition, it converts the coordinates to lat/lon (CRS WGS84 EPSG:4326).

Needed:
    - path_in: path that contains all the ASCII files from WaSiM (temperature and precipitation)
    - file_hgt: DEM, ASCII grid file that contains the altitudes for each grid cell
    - path_out: path where the nc file is going to be saved
    - file_out: name of the nc file that is created, containing both variables (temp and prcp)

"""

import numpy as np
import os
from netCDF4 import Dataset
import datetime
import gdal
import osr
import re

# ===== INPUTS & DEFINITIONS =====
# Path containing the output grids from WaSiM (monthly resolution)
path_in = os.path.abspath(r'E:\DIRT_X\OGGM\meteo_wasim\monthly_values\INCA_clipped_1969-2019')
# Sample path for temperature files from WaSiM (mean monthly values)
pat_temp = re.compile('tempalm_[0-9]{6}.mit')
# Sample path for precipitation files from WaSiM (total monthyl values)
pat_prcp = re.compile('precalm_[0-9]{6}.sum') # sample path for precipitation files

# File with the DEM used in WaSiM (spatial resolution must match meteo grids)
file_hgt = os.path.join(r'E:\DIRT_X\Gepatsch\Gepatschalm\wasim_alm\input', 'dem.asc')
# Path where the converted output files will be saved
path_out = os.path.abspath(r'E:\DIRT_X\Gepatsch\Coupling')
# Name of the file (*.nc) containing the meteorological data to be read by OGGM
file_out = os.path.join(path_out, 'monthly_meteo_inca.nc')
# Define CRS of the input data (used in WaSiM). Here: EPSG:31254 - MGI / Austria GK West – Projected 
epsg_in = 31254
# Define starting date of the files
basedate = datetime.datetime(1969,1,1,0,0,0) 


# First, open the DEM grid with gdal to extract shape and x, y coordinates
# the shape and size is the same for all grids in WaSiM
# https://gis.stackexchange.com/questions/70458/convert-timeseries-stack-of-gtiff-raster-to-single-netcdf    
ds = gdal.Open(file_hgt)
a = ds.ReadAsArray()
nlat, nlon = np.shape(a)
b = ds.GetGeoTransform()

# Set the coordinate reference system (source and destination) and creates 2 arrays with the x and y coordinates
# https://gis.stackexchange.com/questions/177061/ascii-file-with-latitude-longitude-and-data-to-geotiff-using-python
src_crs = osr.SpatialReference()
src_crs.ImportFromEPSG(epsg_in) # source coordinate system
dest_crs = osr.SpatialReference()
dest_crs.ImportFromEPSG(4326) # destination coordinate system
tx = osr.CoordinateTransformation(src_crs, dest_crs) # transform coordinates

cellsize, llx, ury, ncols, nrows = b[1], b[0], b[3], nlon, nlat
lly = ury - cellsize * nrows
(ll_lat, ll_lon, llz) = tx.TransformPoint(lly, llx)
(ur_lat, ur_lon, lrz) = tx.TransformPoint(ury, llx + ncols * cellsize)
mem = gdal.GetDriverByName('MEM')
dest_ds = mem.Create('', nlon, nlat, 1, gdal.GDT_Float32)
dest_geo = (ll_lon, (ur_lon - ll_lon) / ncols, b[2], ur_lat, b[4], (ur_lat - ll_lat) / nrows)
dest_ds.SetGeoTransform(dest_geo)
dest_ds.SetProjection(dest_crs.ExportToWkt())

lon = np.arange(nlon)*dest_geo[1]+dest_geo[0] # lon obtained from x coordinates
lat = np.arange(nlat)*dest_geo[5]+ll_lat # lat obtained from y coordinates

nco = Dataset(file_out ,'w',format='NETCDF4') # create the netCDF file

chunk_lon = nlon
chunk_lat = nlat
chunk_time = 50
    
# Add dimensions
dim_lon = nco.createDimension('lon', nlon)
dim_lat = nco.createDimension('lat', nlat)
dim_time = nco.createDimension('time', None) # in this case the time dimension is "unlimited"
for dim in nco.dimensions.items():
    print(dim)
    
# Add variables
nc_lon = nco.createVariable('lon', 'f4', ('lon'))
nc_lon.units = 'degrees_east'
nc_lon.long_name = 'longitude'
nc_lon[:] = dest_geo[0] + dest_geo[1] / 2 + np.arange(nlon) * dest_geo[1]

nc_lat = nco.createVariable('lat', 'f4', ('lat'))
nc_lat.units = 'degrees_north'
nc_lat.long_name = 'latitude'
nc_lat[:] = np.flipud(ll_lat + dest_geo[5] / 2 + np.arange(nlat) * dest_geo[5])
# lly + cellsize / 2 + np.arange(nrows) * cellsize
nc_time = nco.createVariable('time', 'f4', ('time'))
nc_time.units = 'days since 1969-01-01'
nc_time.calendar = 'proleptic_gregorian'

nc_hgt = nco.createVariable('hgt', 'f4', ('lat', 'lon'))
nc_hgt.units = 'm'
nc_hgt.long_name = 'altitude'
nc_hgt[:] = a

# Create container variable 
crso = nco.createVariable('crs','i4')
crso.long_name = 'Lon/Lat Coords in WGS84'
crso.grid_mapping_name='latitude_longitude'
crso.epsg_code = 'EPSG:4326'
crso.longitude_of_prime_meridian = 0.0
crso.semi_major_axis = 6378137.0
crso.inverse_flattening = 298.257223563

# Create variable for temperature data, with chunking
var_temp = nco.createVariable('temp', 'f4',  ('time', 'lat', 'lon'), 
   zlib=True,chunksizes=[chunk_time,chunk_lat,chunk_lon],fill_value=-9999)
var_temp.units = 'degC'
#tmno.add_offset = 0.00
var_temp.long_name = 'monthly temperature'
var_temp.standard_name = 'air_temperature'
var_temp.grid_mapping = 'crs'
var_temp.set_auto_maskandscale(False)

# create variable for precipitation data, with chunking
var_prcp = nco.createVariable('prcp', 'f4',  ('time', 'lat', 'lon'), 
   zlib=True,chunksizes=[chunk_time,chunk_lat,chunk_lon],fill_value=-9999)
var_prcp.units = 'mm'
#tmno.add_offset = 0.00
var_prcp.long_name = 'monthly precipitation'
var_prcp.standard_name = 'precipitation_amount'
var_prcp.grid_mapping = 'crs'
var_prcp.set_auto_maskandscale(False)

nco.Conventions = 'CF-1.6'

itime_t=0
itime_p=0

#step through data, writing time and data to NetCDF
for root, dirs, files in os.walk(path_in):
    dirs.sort()
    files.sort()
    for f in files:
        if re.match(pat_temp,f): # for temperature
            # read the time values by parsing the filename
            year=int(f[8:12])
            mon=int(f[12:14])
            date=datetime.datetime(year,mon,1,0,0,0)
            dtime=(date-basedate).total_seconds()/86400.
            nc_time[itime_t]=dtime
            temp_path = os.path.join(root,f)
            temp=gdal.Open(temp_path)
            tmn=temp.ReadAsArray()
            var_temp[itime_t,:,:]=tmn
            itime_t=itime_t+1
        elif re.match(pat_prcp,f): # for precipitation
            # read the time values by parsing the filename
            year=int(f[8:12])
            mon=int(f[12:14])
            date=datetime.datetime(year,mon,1,0,0,0)
            dtime=(date-basedate).total_seconds()/86400.
            nc_time[itime_p]=dtime
            prcp_path = os.path.join(root,f)
            prcp=gdal.Open(prcp_path)
            prn=prcp.ReadAsArray() 
            var_prcp[itime_p,:,:]=prn
            itime_p=itime_p+1

nco.close()