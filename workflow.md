## Coupling scheme workflow
This text summarizes the main steps needed to follow for running WaSiM-OGGM
_____

1. Run WaSiM

Daily (or sub-daily) simulations. Monthly mean and sum grids for temperature and precipitation have to be created, respectively. This is done in the interpolation section of the control file.

Variable      | Write code  | Output                
:-----------: |:-----------:| :--------------------:
temperature   | 73          | mean monthly values   
precipitation | 33          | total monthly values 

_____
2. Convert grids to nc files

Use the file: [Convert WaSiM grids to netCDF files](../main/convert_raster_to_nc.py) \
Results: one netCDF file (i.e.: 'monthly_meteo.nc') containing monthly values of temperature and precipitation.

_____
3. Adjust monthly meteo data for OGGM requirements

In case the glaciers do not belong to the reference glaciers, an adjustement of the meteorological input data for OGGM is required. This is needed because the parameters are obtained through interpolation from neighboring glaciers, which were already calibrated using default climate datasets.
The following steps are required:

* Run OGGM with default climate dataset (CRU)
* Run OGGM with monthly values of temperature and precipitation prepared during step 2.
* Run the script [Adjust meteorological data for OGGM](../main/adjust_meteo.py) 
* The user dataset (x_user) is adjusted (x_adj) to the default dataset (y_def)
  
  * Linear regression for temperature: slope and interception are obtained
  * Scaling factor for precipitation: user dataset is adjusted with a scaling factor
* In both cases, the adjustement is based on the mean monthly values of the variables.
* Modify the netCDF file created in 2. to adjust the values according to x_adj with [CDO (Climate Data Operator)](https://code.mpimet.mpg.de/projects/cdo), running in Linux:

  * cdo -expr, 'temp = (temp + b)/a' inputnc4 output_temp  <!--- linear regression for temperature -->
  * cdo -expr, 'prcp = prcp * alfa' inputnc4 output_prcp   <!--- scaling factor for precipitation -->
  * cdo -expr, 'hgt = hgt' inputnc4 output_hgt             <!--- create a nc file with the reference altitudes (needed for merging all files then -->
  * cdo merge output_temp output_prcp output_hgt output_monthly_data.nc <!--- merge the outputs with the adjusted meteo data and the height and create a new nc file which will be used finally by OGGM -->
_____
4. Run OGGM with adjusted meterological data from WaSiM


_____
5. Process annual glacier outlines


_____
6. Re-run WaSiM with initial states from OGGM outputs


_____
7. Optimize WaSiM results (calibration)
