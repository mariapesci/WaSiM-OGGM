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
Use the file: [Convert WaSiM grids to netCDF files](../blob/master/convert_)

_____
3. Adjust monthly meteo data for OGGM requirements


_____
4. Run OGGM


_____
5. Process annual glacier outlines


_____
6. Re-run WaSiM with initial states from OGGM outputs


_____
7. Optimize WaSiM results (calibration)
