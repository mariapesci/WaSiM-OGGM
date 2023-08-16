# -*- coding: utf-8 -*-
"""
@author: María Pesci

This script reads the glacier areas and thicknesses created after running OGGM and merged them together in one file
containing the ice thickness per year.

Needed:
    - path_thick: path containing the annual (per glacier) distributed thickness created from OGGM and the time series
      with the total glacier volume from OGGM ('vol.csv')
    - path_out: path in which the merged outputs for running WaSiM (containing glacier thicknesses) will be saved
    - dem_wasim: DEM used in the WaSiM model

"""

import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import os
import lmfit
import scipy as sp
from scipy.optimize import minimize_scalar
import copy
from pandas.plotting import table
import glob
import geopandas as gpd
import rioxarray
from rioxarray.merge import merge_arrays
import xarray
import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.warp import calculate_default_transform, reproject, Resampling

import wasimlib

path_thick = r'E:\DIRT_X\Gepatsch\Coupling\ice_thickness_vol\oggm_v1.6\thick_dyn_spinup_calibrated\from_oggm'
path_out = r'E:\DIRT_X\Gepatsch\Coupling\ice_thickness_vol\oggm_v1.6\thick_dyn_spinup_calibrated\for_wasim' 

# 1. Read the glacier thickness per glacier and reproject/resample it to the project's CRS (https://rasterio.readthedocs.io/en/latest/topics/reproject.html)
thick_tif = glob.glob(path_thick + '\*.tif')

crs_out = 'EPSG:31254'
dem_wasim = r'E:\DIRT_X\Gepatsch\Gepatschalm\GIS_Gepatschalm\dem_wasim_100m.tif'

with rasterio.open(dem_wasim) as src:
    transform, width, height = calculate_default_transform(
        src.crs, crs_out, src.width, src.height, *src.bounds)
    kwargs = src.meta.copy()
    kwargs.update({
        'crs' : crs_out,
        'transform': transform,
        'width': width,
        'height': height,
        'nodata': -9999.})
    
for tf in thick_tif:
        with rasterio.open(os.path.join(path_thick, tf)) as glc_thick:
            transform2, width2, height2 = calculate_default_transform(
                glc_thick.crs, crs_out, glc_thick.width, glc_thick.height, *glc_thick.bounds)
            #print(glc_thick.nodata)
            with rasterio.open(os.path.join(path_out, tf.split('glac_thick_')[1].split('.tif')[0]+ '_rep.tif'), 'w', **kwargs) as dst:
                #print(dst.nodata)
                dst.nodata = -9999.                
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(glc_thick,i),
                        destination=rasterio.band(dst,i),
                        src_transform=glc_thick.transform,
                        src_crs=glc_thick.crs,
                        dst_transform=src.transform,
                        dst_crs = crs_out,
                        nodata = -9999.,
                        resampling=Resampling.nearest)

#%%
# 2. Merge all glacier files together (https://corteva.github.io/rioxarray/stable/examples/merge.html)
# Define initial year
init_year = 1969
# And period (in years) for which the files will be created
years = np.arange(init_year+1, init_year+51) 

for year in years:
    thick_merged = glob.glob(path_out + '\*'+str(year)+'_rep.tif')
    list_for_merging = []
    for i in range(len(thick_merged)):
        file = rioxarray.open_rasterio(thick_merged[i])[0]                                       
        list_for_merging.append(file)

    merged = merge_arrays(dataarrays=list_for_merging)    
    merged.rio.to_raster(os.path.join(path_out, 'merged_'+str(year)+'.tif'))

   
#%%
# 3. Get the "volume" grid from area and thickness
vol_oggm = pd.read_csv(os.path.join(path_thick, 'vol.csv'), index_col='time')
vol_index = vol_oggm.index.astype(int)
vol_oggm.set_index(vol_index, inplace=True)

for year in years:
    glac_cells = rioxarray.open_rasterio(os.path.join(path_thick, 'glaciercells'+str(year)+'.asc'))[0].values
    glac_thick = rioxarray.open_rasterio(os.path.join(path_out, 'merged_'+str(year)+'.tif'))[0].values

    min_thick = np.min(glac_thick[np.nonzero(glac_thick)])
    vol = np.nan_to_num(glac_cells, copy=True)
    vol_conserv = vol.copy()
    sum_vol = []
    sum_vol_conserv = []
    vol_i = 0
    for r in range(glac_cells.shape[0]):
        for c in range(glac_cells.shape[1]):
            #print(r,c)
            if glac_cells[r,c] != -9999:
                if glac_thick[r,c] > 0:
                    vol[r,c] = glac_cells[r,c] * glac_thick[r,c]
        
            if vol[r,c] != -9999 and vol[r,c] < min_thick:
                vol[r,c] = min_thick
            if vol[r,c] != -9999:
                vol_i = vol[r,c] * (100 * 100)
                vol[r,c] = vol[r,c] * (100 * 100)
                sum_vol.append(vol_i)
        

    # Conserve volume (from OGGM)      
    vol_oggm_year = vol_oggm.loc[year].vol_m3 
    print(year, vol_oggm_year, sum(sum_vol))
    scal_vol = vol_oggm_year/sum(sum_vol)
    #print(scal_vol)
    for r in range(vol.shape[0]):
        for c in range(vol.shape[1]):
            if vol_conserv[r,c] != -9999:
                vol_conserv[r,c] = scal_vol * vol[r,c] / (100 * 100)


    ncols = vol_conserv.shape[1]
    nrows = vol_conserv.shape[0]
    # xll and yll corner have to be the same as used in WaSiM grids
    xllcorner = 22597.5 
    yllcorner = 185062.5 
    cellsize = 100.000000000000
    nodata_value = -9999

    # Finally, a grid containing the ice thickness for the entire catchment is created per year
    vol_init = os.path.join(path_out, 'glaciercells'+str(year)+'.asc')
    write_head = "ncols\t%d\nnrows\t%d\nxllcorner\t%f\nyllcorner\t%f\ncellsize\t%f\nnodata_value\t%f" % \
                (ncols,nrows,xllcorner,yllcorner,cellsize,nodata_value)
    np.savetxt(vol_init,vol_conserv,fmt="%f",header=write_head,comments="")
