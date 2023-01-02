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
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import copy
from netCDF4 import Dataset, date2num
import geopandas as gpd
from geocube.api.core import make_geocube
import rioxarray
import salem
import shapely.geometry as shpg
from shapely.ops import unary_union
from oggm import utils
from oggm.core.flowline import FileModel
import fiona
import rasterio
from rasterio.transform import from_origin


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


           
def from_width_to_polygon(gdir, year, crs_out, path_for_shp, save_shp=True):
    """
    This function delineates the polygon (glacier outline) from the catchment
    widths, along the different flowlines. It also corrects the total area, to 
    make it comparable to the calculated one by OGGM. Finally, it saves the
    polygon (glacier outline) per glacier and year into a shp file, with the
    specified CRS.

    Parameters
    ----------
    gdir : GlacierDirectory
        the glacier directory
    year : int
        year for which the outline should be created
    crs_out : str
        coordinates reference system. Example: 'EPSG:31254'
    path_for_shp : str
        path in which the glacier outline (shapefile) will be saved
    save_shp : bool, optional
        if the shapefile should be saved. The default is True.

    Returns
    -------
    the glacier outline for the specified year.

    """
    list_pol = []
    areas_gl = []
    smap = salem.Map(gdir.grid, countries=False, nx=gdir.grid.nx)
    with utils.ncDataset(gdir.get_filepath('gridded_data')) as nc:
        topo = nc.variables['topo'][:]
    smap.set_topography(topo)
    model = FileModel(gdir.get_filepath('model_geometry'))
    model.run_until(year)
    crs = gdir.grid.center_grid
    clss = model.fls[::-1]
    for i, fl in zip(range(len(clss)), clss):
        a_cur, b_cur, a_b = [], [], []
        area_fl = fl.area_m2
        for wi, cur, (n1, n2) in zip(clss[i].widths, clss[i].line.coords, 
                         clss[i].normals):

            _l = shpg.LineString([shpg.Point(cur+wi/2.*n1), shpg.Point(cur+wi/2.*n2)])
            aiw = shpg.Point(cur+wi/2.*n1)
            biw = shpg.Point(cur+wi/2.*n2)
            ai = salem.transform_geometry(aiw, crs=crs, to_crs=crs_out)
            bi = salem.transform_geometry(biw, crs=crs, to_crs=crs_out)
            point_a = (ai.x, ai.y)
            point_b = (bi.x, bi.y)
            a_cur.append(point_a)
            b_cur.append(point_b)
            a_b = a_cur+ b_cur[::-1]
            smap.set_geometry(_l, crs=crs, color='r', linewidth=0.6, zorder=50)
            smap.set_geometry(fl.line, crs=crs, color='r', linewidth=0.6, zorder=50)
        pol = shpg.Polygon(a_b) # A) Option to construct normal polygon from the points
        list_pol.append(pol.buffer(0))
        areas_gl.append(area_fl)
    multi_pol = unary_union(list_pol) 

    # Conserve area
    buffer_area = (sum(areas_gl) - multi_pol.area) / multi_pol.length
    glac_pol = multi_pol.buffer(buffer_area)
    glac_pol.area
    smap.set_geometry(glac_pol, crs=crs_out, linewidth=0.6, zorder=50)

    # Save polygon as shapefile:
    if save_shp is True:
        shpef = path_for_shp + gdir.rgi_id + '_' + str(year) + '.shp'
        schema =  {'geometry': 'Polygon', 'properties': {'id': 'int'},}
        with fiona.open(shpef, 'w', 'ESRI Shapefile', schema, crs=crs_out) as output:
            output.write({
                'geometry': shpg.mapping(glac_pol),
                'properties': {'id': 'gep'},
                })

    return

# Rasterize the shapefile with the fraction of glacierized cells and glacier codes
def rasterize_shp(path_shp, path_raster_dem, shapefile_name, year):
    """
    This function reads a shp file and converts it into raster (ASCII format)
    to be read by WaSiM. The shapefiles contain the grid cells of the model and are:
        - glaciercells: grids with fraction of glacierized cells
        - glaciercodes: grids with the corresponding code of each glacier (it always
          has to match the glaciercells and the subbasins codes)

    Parameters
    ----------
    path_shp : str
        path with the location where the shapefiles are
    path_raster_dem : str
        path with the location of the DEM used by WaSiM (containing grid cell information)
    shapefile_name : str
        either glaciercells (for fraction of glacierized cells) or glaciercodes (for
        codes of each glacierized cell)
    year : int
        year for which the rasters will be created

    """
    raster_create = os.path.join(path_shp, shapefile_name + str(year) + '.asc')
    source_ds = gpd.read_file(os.path.join(path_shp, shapefile_name + str(year) + '.shp'))
    raster = rioxarray.open_rasterio(path_raster_dem, masked=True)
    cube = make_geocube(vector_data=source_ds, like=raster, fill=-9999)
    if shapefile_name == 'glaciercells':
        cube['frac_glac'].rio.to_raster(raster_create)
    elif shapefile_name == 'glaciercodes':
        cube['code'].rio.to_raster(raster_create)
    else:
        print('The name of the shapefile should be '
              'glaciercells for grid containing fraction of glacierized cells'
              'or glaciercodes for grid containing glacier codes.')
    return


def prepare_grids_for_wasim(gdirs, year, crs_out, path_oggm_outlines, path_input_wasim, 
                            dem_wasim, grid_nodes_wasim, subbasins, resol_wasim):
    """
    This function reads the shapefile (glacier outline) per glacier and year and does the
    required operations to convert it into an input for wasim:
        - merge all glaciers together
        - clip the glaciers to the wasim grid
        - prepare the subbasins from wasim so the glaciers match the codes of them
        - get the cells that are contained in the outline and create glaciercells
        and glaciercodes grids
        - rasterize the previous grids according to the dem used in wasim

    Parameters
    ----------
    gdirs : list
        list with the glacier directories
    year : int
        year for which the grids should be created
    crs_out : str
        coordinates reference system. Example: 'EPSG:31254'
    path_oggm_outlines : str
        path to the shapefiles (glacier outlines) created for each glacier and year
    path_input_wasim : str
        path to which the new files will be saved
    dem_wasim : str
        path to the DEM from WaSiM
    grid_nodes_wasim : str
        path to the grid shapefile from WaSiM
    subbasins : str
        path to the subbasins from WaSiM
    resol_wasim : int
        spatial resolution used in WaSiM. It should match the resolution of the DEM

    Returns
    -------
    Annual grids in ASCII format of glaciercells and glaciercodes, ready to be read by WaSiM.

    """
    
    shp_files = glob.glob(os.path.abspath(path_oggm_outlines + '/*'+str(year)+'.shp')) # read outlines from oggm
    # First step: merge all shp of all glaciers into one (concatenate): 
    glac_out = gpd.GeoDataFrame(pd.concat([gpd.read_file(shp) for shp in shp_files]))
    glac_out['id'] = range(len(gdirs))
    # check if there are any empty polygons (geometry = 0)
    for geom in range(len(glac_out)):
        if glac_out['geometry'].iloc[geom] is None:
            glac_out['geometry'].iloc[geom] = shpg.Polygon([])
    
    glac_out.to_file(os.path.join(path_input_wasim, 'glac_out_oggm_'+str(year)+'.shp'))
    
    # Intersect outline of glacier with grid WaSiM
    grid_wasim = gpd.read_file(grid_nodes_wasim) # read wasim grid
    
    glac_clip = gpd.clip(grid_wasim, glac_out) # clip both files
    # Add attribute field to get the fraction of glacierized cell
    glac_clip['glac_cells'] = glac_clip.area / (resol_wasim)
    glac_clip.to_file(os.path.join(path_input_wasim,'glac_cells_'+str(year)+'.shp')) # save glacierized fraction shp
    
    # Join attributes by field value    
    frac_glac = grid_wasim.merge(glac_clip, on='node')
    # Convert pandas DataFrame to GeoDataFrame
    frac_glac_cells = gpd.GeoDataFrame(frac_glac, crs=crs_out, geometry=frac_glac['geometry_x'])


    # Read the file containing the subbasins (same as glacier codes)
    subbasins_wasim = gpd.read_file(subbasins)
    frac_glac_cells['frac_glac'] = -9999
    for si in subbasins_wasim.index:
        for fi in frac_glac_cells.index:
            if subbasins_wasim['geometry'].loc[si].contains(frac_glac_cells['geometry'].loc[fi]):
                frac_glac_cells['frac_glac'].loc[fi] = frac_glac_cells['glac_cells'].loc[fi]

    # Gepatsch    
    frac_glac_cells.drop(columns=['geometry_x', 'geometry_y','glac_cells'], inplace=True)
    # Rofenache
    #frac_glac_cells.drop(columns=['id_x','geometry_x','id_y','geometry_y','glac_cells'], inplace=True)
    # Save new GeoDataFrame with fraction of glacierized cells as shapefile
    frac_glac_cells.to_file(os.path.join(path_input_wasim, 'glaciercells'+str(year)+'.shp'))

    # Now we have to create the file containing the glacier codes
    glac_codes = pd.DataFrame.copy(frac_glac_cells) # make a copy of the shp with the glacier's information
    glac_codes['code'] = -9999 # create a column which will contain the code of the glacier (subbasin id)
    # Iterate through each glacier cell to find to which subbasin id belongs
    for si in subbasins_wasim.index:
        for gi in glac_codes.index:
            if subbasins_wasim['geometry'].loc[si].contains(glac_codes['geometry'].loc[gi]):
                glac_codes['code'].loc[gi] = subbasins_wasim['ID'].loc[si]

    glac_codes.to_file(os.path.join(path_input_wasim, 'glaciercodes'+str(year)+'.shp'))

    rasterize_shp(path_input_wasim, dem_wasim, 'glaciercodes', year)
    rasterize_shp(path_input_wasim, dem_wasim, 'glaciercells', year)
    
    return



def topo_and_mask_per_glacier(gdir, year, path_shp, save_mask=True):
    """
    This function reads the topography from each glacier directory created by OGGM
    (which was saved as 'gridded_data') and the shapefile containing the glacier
    outline created by the function "from_width_to_polygon()" and:
        - saves the topography into a raster 
        - clips the glacier outline (shapefile) with the corresponding topography (raster)
        - saves the glacier masks if specified
        - returns the topo_smoothed and glacier_mask (as in the function given by OGGM:
          "distribute_thickness_per_altitude()") as arrays.

    Parameters
    ----------
    gdirs : list
        list with the glacier directories
    year : int
        year for which the grids should be created
    path_shp : str
        path in which the glacier outline (shapefile) is saved
    save_mask : bool, optional
        if the glacier_mask raster is to be saved. The default is True.

    Returns
    -------
    topo_smoothed : numpy.ndarray
        array with topo_smoothed
    glacier_mask : numpy.ndarray
        array with glacier masks (1=glacier, 0=no glacier)

    """
    
    t = xr.open_dataset(gdir.get_filepath('gridded_data'))
    # Open the gridded data used by OGGM (topography) for getting the extent and resolution
    with utils.ncDataset(gdir.get_filepath('gridded_data')) as nc:
        topo = nc.variables['topo_smoothed'][:]
    # Read the glacier outline for specified year, created previously
    glac_outline = gpd.read_file(path_shp + gdir.rgi_id + '_' + str(year) + '.shp')
    glac_outline = glac_outline.to_crs(t.pyproj_srs)
    
    # Convert numpy array with topo_smoothed data to a raster
    # Define the data extent (min. lon, min. lat, max. lon, max. lat)
    extent = [t.x.values[0], t.y.values[0], t.x.values[-1], t.y.values[-1]]   #t.glacier_ext
    dx = t.x.values[1] - t.x.values[0]
    dy = t.y.values[0] - t.y.values[1]
    transform = from_origin(extent[0], extent[1], dx, dy)
    topo_raster = rasterio.open(os.path.join(path_shp,'topo_raster_'+gdir.rgi_id+'_'+str(year)+'.tif'),
                                'w', driver='GTiff',
                                height=topo.shape[0], width=topo.shape[1],
                                count=1, dtype=str(topo.dtype),
                                crs=t.pyproj_srs,
                                transform=transform)
    topo_raster.write(topo,1)
    topo_raster.close()
    
    # Intersect/clip the raster with topo_smoothed data and the shp containing the outline
    #topo_raster_open = rioxarray.open_rasterio(os.path.join(path_for_shp, topo_smoothed_raster))

    # clip the glacier outline to the raster extent
    with rasterio.open(os.path.join(path_shp, 'topo_raster_'+gdir.rgi_id+'_'+str(year)+'.tif')) as src:
        out_image, out_transform = rasterio.mask.mask(src, glac_outline['geometry'], crop=False)
        out_meta = src.meta
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})
    if save_mask:
        with rasterio.open(os.path.join(path_shp, 'glac_masked_'+gdir.rgi_id+'_'+str(year)+'.tif'), "w", **out_meta) as dest:
            dest.write(out_image)   

    topo_smoothed = out_image[0].copy()
    glacier_mask = out_image[0].copy()
    glacier_mask[glacier_mask != 0] = 1
    
    crs_mask = t.pyproj_srs
    
    return topo_smoothed, glacier_mask, crs_mask


def glac_thickness(gdir, year, path_out, glacier_mask):
    t = xr.open_dataset(gdir.get_filepath('gridded_data'))
    # Open the gridded data used by OGGM (topography) for getting the extent and resolution
    with utils.ncDataset(gdir.get_filepath('gridded_data')) as nc:
        topo = nc.variables['topo_smoothed'][:]
            
    cl_s = gdir.read_pickle('inversion_output')
    fl_s = gdir.read_pickle('inversion_flowlines')
    hs, ts, vs, ws, xs, ys = [],[],[],[],[], []
    for cl, fl in zip(cl_s, fl_s):
        hs = np.append(hs, fl.surface_h)
        ts = np.append(ts, cl['thick'])
        vs = np.append(vs, cl['volume'])
        ws = np.append(ws, fl.widths_m)
        
        try:
            x, y = fl.line.xy
        except AttributeError:
            x = fl.surface_h * 0 - 1
            y = fl.surface_h * 0 - 1
            
        xs = np.append(xs, x)
        ys = np.append(ys, y)

        
    init_vol = np.sum(vs)

    thick = glacier_mask + np.NaN
    for y in range(thick.shape[0]):
        for x in range(thick.shape[1]):
            phgt = topo[y,x]
            starth = 100.
            while True:
                starth +=10
                pok = np.nonzero(np.abs(phgt - hs) <= starth)[0]
                if len(pok) != 0:
                    break
            sqr = np.sqrt((xs[pok]-x)**2 + (ys[pok]-y)**2)
            pzero = np.where(sqr==0)
            if len(pzero[0]) == 0:
                thick[y,x] = np.average(ts[pok], weights=1/sqr)
            elif len(pzero[0]) == 1:
                thick[y,x] = ts[pzero]
            else:
                raise RuntimeError('wrong')

    dx = gdir.grid.dx
    dy = gdir.grid.dy
    # Re-mask
    utils.clip_min(thick, 0, out=thick)
    thick[glacier_mask == 0] = np.NaN
    assert np.all(np.isfinite(thick[glacier_mask == 1]))
    # Conserve volume
    tmp_vol = np.nansum(thick * dx**2)
    thick *= init_vol / tmp_vol

    thick[glacier_mask == 0] = -9999. # make nodata as -9999.0
    extent = [t.x.values[0], t.y.values[0], t.x.values[-1], t.y.values[-1]]   #t.glacier_ext
        
    glac_thick = rasterio.open(os.path.join(path_out,'glac_thick_'+gdir.rgi_id+'_'+str(year)+'.tif'),
                               'w', driver='GTiff',
                               height=thick.shape[0], width=thick.shape[1],
                               count=1, dtype=str(thick.dtype),
                               crs=t.pyproj_srs,
                               nodata= -9999.,
                               transform=from_origin(extent[0], extent[1], dx, -dy))
    glac_thick.write(thick,1)
    glac_thick.close()
    
  
    return thick

#def resample_raster(gdir, year, path_raster, path_glaccells, wasim_dx):
    # Reproject to Wasim's CRS and pixel resolution




    # scale image transform
    # transform = thick_raster.transform * thick_raster.transform.scale(
    #     (thick_raster.width / data.shape[-1]),
    #     (thick_raster.height / data.shape[-2]))


# def read_init_results(gdir, path_init_results, init_year):
#     """
#     This function reads the results (tar file) of the initialization run (Eis et al., 2021),
#     extract the files, unpickle the results of the initialization year for each glacier,
#     selects the best fitness value and the corresponding model geometry.

#     Parameters
#     ----------
#     gdirs : GlacierDirectory
#         the glacier directory (not as a list)
#     path_init_results : str
#         path with the location of the results of the initialization run 
#     init_year: int
#         starting year of the initialization run, which is also the starting year
#         of the actual run
        
#     Returns
#     -------
#     file : str
#         name of the results file with the best fit

#     """
#     path_init = os.path.join(path_init_results, 'per_glacier')

#     region = gdir.rgi_id.split('.')[0]
#     fp = path_init+'/'+region+'/'+region+'.00'+'/'+gdir.rgi_id+'.tar.gz'
#     # Define the path where the "tar" files of the initialization run are stored
#     init_dir = path_init_results+'/'
#     # Extract the "tar" files
#     with tarfile.open(fp, 'r') as tf:
#         tf.extractall(os.path.dirname(init_dir))
#     # Read the results folder (pickle) for each of the glaciers. It starts with the init_year
#     result_init = init_dir+gdir.rgi_id+'/result'+str(init_year)+'.pkl'
#     # Unpickle the results folder
#     with gzip.open(result_init, 'rb') as f:
#         result = pickle.load(f)
#     # Select the glacier geometry with the smallest fitness value, which will be the 
#     # initial geometry for the initialization in init_year
#     best_fitness = result.sort_values('fitness', ascending=True)['fitness'].iloc[0]
#     temp_bias_init = result.sort_values('fitness', ascending=True)['temp_bias'].iloc[0]
#     time_init = result.sort_values('fitness', ascending=True)['time'].iloc[0]
#     # File with the best fit, that will be used afterwards
#     filename_best_fit = str(init_year)+'_past_*'+str(temp_bias_init)+'_'+str(time_init)+'.nc'
#     file = glob.glob(os.path.join(init_dir, gdir.rgi_id, str(init_year), 'model_geometry'+filename_best_fit))[0]
#     shutil.copy(file, gdir.dir)
    
#     return file

