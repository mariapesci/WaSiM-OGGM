"""
@author: María Pesci

This script contains the settings and functions for performing the automatic calibration
of the third (WaSiM) run of the coupling scheme. It does not produce any output, it has
to be run with 'spotpy_calibration.py'.

An example of the control file required for the coupling calibration can be found in 'ctrl_coupling_master.txt'
"""

import sys
import os
import pandas as pd
import numpy as np
from spotpy.parameter import Uniform
from spotpy.objectivefunctions import rmse
import datetime
import shutil
import glob

import wasimlib
import funct_coupling

def update_ctrl_file(ifile, ofile, variables_str, variables_val):
    """Updates a master control file through replacing variables by values
        specified by the user
        
        Parameters
        ----
        ifile : path to input file (master file)
        ofile : path to output file (variables will be replaced by values)
        variables_str : names of the variables to be replaced
        variables_val : the parameter set that is exported to the control file
                
    """       

    if len(variables_str) != len(variables_str) or len(variables_str) == 0 or len(variables_val) == 0:
        print('Invalid parameters!')
        return False
        # prepare string templates
    list_var_ids = list()
    list_var_val = list()
    list_var_upd = list()
    n_v = len(variables_str)
    for i in range(0,n_v):
        cmpstr = '§%s§' % variables_str[i]
        list_var_ids.append(cmpstr)
        valstr = '%.4f' % variables_val[i]
        list_var_val.append(valstr)
        list_var_upd.append(False)

    n_r = 0        
    try:
        fid = open(ofile,"w+", encoding="utf-8")
        with open(ifile, encoding="utf-8", errors='ignore') as f:
            for line in f:
                line_updated = str(line)
                for i in range(0,len(variables_str)):               
                    line_updated = line_updated.replace(list_var_ids[i], list_var_val[i])
                if line_updated != line:
                    n_r += 1
                # lookup multiple operators
                split1 = line_updated.split('<<')
                if len(split1) == 2:
                    if len(split1[0]) < len(line):
                        # in this case, multiple operators might have been detected
                        split2 = split1[1]
                        if len(split2) > 0:
                            # get operator
                            operator = split2[0]
                            #print(operator)
                            split3 = split2[1:-1].split('>>')
                            if len(split3) == 2:
                                if operator == '*':
                                    # multiplication of all subsequent numbers
                                    try:
                                        argument = float(split3[0])
                                        #check for semicolons
                                        semi = split3[1].split(';')
                                        line_of_numbers = semi[0]
                                        add_semi = False
                                        if len(semi) > 1:
                                            add_semi = True
                                        split4 = line_of_numbers.split(' ')
                                        number_list = list()
                                        for substr in split4:
                                            try:
                                                number = float(substr)
                                                number_list.append(number * argument)
                                            except:
                                                pass
                                        if len(number_list) > 0:
                                            # in this case we have valid data
                                            # new line update
                                            line_updated = split1[0]
                                            for ni in number_list:
                                                line_updated = line_updated + '  ' + str(ni)
                                            if add_semi:
                                                line_updated = line_updated + ' ;'
                                            line_updated = line_updated + '\n'
                                    except:
                                        pass

                                
                fid.write(line_updated)
        fid.close()
    except:
        print('An error occured while processing the files')
        return False

    return True
    
class wasim_coupling_model(object):
    def __init__(self, coupling_settings, ctrl_master, ctrl_cal, 
        param_var, param_val, sim_runoff, sim_glmb, target_id_runoff=None, target_id_glac=None):
        """ This functions sets all relevant attributes of the wasim object.
        
        Parameters
        ----
        coupling_settings : class object
                            class containing all input data and variables needed for running wasim and the calibration
        ctrl_master       : str
                            location of the control file with parameters defined as variables (§ and &)
        ctrl_cal          : str
                            location of the control file that will be generated in each calibration run (for entire sim period)
        param_var         : list str
                            list with the name of the parameters that will change during calibration
        param_val         : list num
                            list with the value of the parameters that changes in each calibration step
                            (param_var and param_val are defined within spotpy_calibration.py/simulation)
        sim_runoff        : str
                            name of the file containing simulated runoff from WaSiM (e.g. 'qgkoalm')
        sim_glmb          : str
                            name of the file containing simulated mass balance of the glaciers from WaSiM (e.g. 'glmb2alm')
        target_id_runoff  : str
                            number of the subcatchment for which observed and simulated runoff is compared (calibrated)
        target_id_glacier : str
                            number or name of the subcatchment for which observed and simulated glacier mass balance is compared (calibrated)

        """
        self.config = coupling_settings
        self.ctrl_master = ctrl_master
        self.ctrl_cal = ctrl_cal
        self.sim_runoff = sim_runoff
        self.sim_glmb = sim_glmb
        self.param_var = param_var
        self.param_val = param_val
        self.target_id_runoff = target_id_runoff
        self.target_id_glac = target_id_glac

        #print('Copy initial states from %s to %s...' % (self.init_file_in, self.init_file_out))

    def run(self):
        """Calls the WaSiM program and passes appropriate arguments.
        
        Parameters
        ----
        """
        wor_dir          = self.config.working_dir
        path_oggm        = '%s%s' % (wor_dir, self.config.init_glac_dir)
        path_input       = '%s%s' % (wor_dir, self.config.input_dir)
        path_init_states = '%s%s' % (wor_dir, self.config.init_states_dir)
        path_meteo       = '%s%s' % (wor_dir, self.config.meteo_dir)
        ctrl_new         = '%s%s' % (wor_dir, self.config.ctrl_new)
        path_out         = '%s%s' % (wor_dir, self.config.output_dir)
        path_results     = '%s%s' % (wor_dir, self.config.results_dir)

        list_steps = np.arange(self.config.year1, self.config.year2, self.config.step)
        delete_outputs = True
        save_files = True

        # 1° create a copy of the control file with variables defined by § for parameters to be calibrated
        # and & for parameters to be updated annually during the coupling 
        # ctrl_master copy to ctrl_cal:
        os.chdir(wor_dir)
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        shutil.copyfile(self.ctrl_master, self.ctrl_cal)
        update_ctrl_file(self.ctrl_master, self.ctrl_cal, self.param_var, self.param_val)

        # 2° run WaSiM annually, updating glacier grids:
        # at the beginning, readgrids = 0 means that no storage grids are read, since the model is initialized now
        # from 2nd year onwards, the outputs of the first year are saved as initial states for the following year
        # and readgrids = 1, so WaSiM knows that it needs to read storage grids (init. files)
        for ii, yi in enumerate(list_steps):
            print('%s' % yi)
            if yi == self.config.year1:
                readgrids = str(0)
            else:
                readgrids = str(1)
            startyear = str(yi)    
            if ii == len(list_steps)-1:
                endyear_i = list_steps[ii] + 1
            else:
                endyear_i = list_steps[ii+1]
            endyear = str(endyear_i)
            
            init_directory = path_init_states + 'init_states_' + startyear + '/'
            if not os.path.exists(init_directory):
                os.makedirs(init_directory)
            out_directory = wor_dir + 'output_' + endyear + '/'
            if not os.path.exists(out_directory):
                os.makedirs(out_directory)
            
            # Copy the glacier grids (from OGGM output) into the initial states directory (for the corresponding year)
            glc_i = self.config.glaccells_grd % (endyear_i) # change names of grd files to match startyear!!!
            glid_i = self.config.glaccodes_grd % (endyear_i)
            shutil.copyfile('%s%s' % (path_oggm, glc_i), '%s%s' % (init_directory, glc_i))
            shutil.copyfile('%s%s' % (path_oggm, glid_i), '%s%s' % (init_directory, glid_i))
            # And re-start the mass balance grids
            glmb_i = self.config.glmb_grd
            glmb_old_i = self.config.glmbold_grd
            if yi > self.config.year1:
                shutil.copyfile('%s%s' % (path_oggm, glmb_i), '%s%s' % (init_directory, glmb_i))
                shutil.copyfile('%s%s' % (path_oggm, glmb_old_i), '%s%s' % (init_directory, glmb_old_i))
                
           
            b_i = self.config.va_sca_factor
            f_i = self.config.va_exp_factor
                
            list_var = ['STARTYEAR', 'ENDYEAR', 'READGRIDS', 'B',  'F',   
                        'INIT_DIRECTORY', 'OUTPUT_DIRECTORY', 'INPUT_DIRECTORY', 'METEO_DIRECTORY']
            list_val = [ startyear,   endyear,   readgrids,   b_i,  f_i,   
                        init_directory,   out_directory,      path_input,        path_meteo      ]
            list_tp  = [ False,       False,     False,       True, True,  
                        False,            False,              False,             False           ] # if False, the variable is a string
            shutil.copyfile(self.ctrl_cal, ctrl_new)
            wasimlib.update_ctrl_file(self.ctrl_cal, ctrl_new, list_var, list_val, list_tp)
            
            #current_path = os.getcwd()
            os.chdir(wor_dir)
            
            # Call wasim to run annually
            wasimlib.call_wasim(['%s' % wor_dir + 'wasim/wasimvzo64.exe', '%s' % ctrl_new], verbose=False)
            # Creates new folder for initialization of the next year, where the current output will be moved
            next_init_directory = path_init_states + 'init_states_' + endyear + '/'            
            if not os.path.exists(next_init_directory):
                    os.makedirs(next_init_directory)
            
            src_files = os.listdir(out_directory)
            # Copy files from current output directory to following initialization directory
            for file_name in src_files:
                full_file_name = os.path.join(out_directory, file_name)
                shutil.copy(full_file_name, next_init_directory)
                # and to the final output with the year of generation
                # useful for concatenate all files at the end of the simulation period
                shutil.copy(full_file_name, path_out)
                os.rename(os.path.join(path_out, file_name), os.path.join(path_out, endyear + '_' + file_name))
                
            # Delete annual outputs if required
            if delete_outputs:
                shutil.rmtree(out_directory)

        # Concatenate all results
        list_files = glob.glob(path_out+'\*.b'+str(self.config.year))
        list_names = []
        for file in list_files:
            if endyear in file:
                file_name = file.split('\\')[1].split(endyear+'_')[1]
                list_names.append(file_name)

        df = pd.DataFrame()
        dict_name = {name: df for name in list_names}

        for key, value in dict_name.items():
            for yi in list_steps:
                df_yi = wasimlib.read_time_series(glob.glob(os.path.join(path_out, str(yi+1) +'_'+key))[0], drop_stats=True)[0]
                value = value.append(df_yi)
            dict_name[key] = value
            if save_files:
                value.to_csv(os.path.join(path_results, key.split('.b')[0] + '.csv'))
        
        shutil.rmtree(path_out)
        
        
        runoff = pd.read_csv(os.path.join(path_results, self.sim_runoff+'.csv'), index_col=0, parse_dates=True)[self.target_id_runoff]
        glmb = pd.read_csv(os.path.join(path_results, self.sim_glmb+'.csv'), index_col=0, parse_dates=True)[self.target_id_glac]
        glmb_annual = glmb.resample('AS-SEP', origin='end').last() #.resample('AS-OCT', origin='end').last().diff()
        
        return runoff, glmb_annual

class coupling_settings(object):
    catchment              = 'Gepatschalm' 

    year                   = '70'
    year1                  = 1969 
    year2                  = 1999 
    step                   = 1
    va_sca_factor          = 1.0
    va_exp_factor          = 1.0

    init_glac_dir          = 'init_states_glac_vol_oggm1.6/dyn_spinup/'
    input_dir              = 'input/'
    init_states_dir        = 'initial_states/'
    meteo_dir              = 'meteo_INCA_kNN/'
    ctrl_master            = 'wasim/ctrl_coupling_master.txt'
    ctrl_cal               = 'wasim/ctrl_coupling_cal.txt'
    ctrl_new               = 'wasim/ctrl_coupling_new.txt'
    output_dir             = 'output/'
    results_dir            = 'results_coupling/'

    working_dir        = r'E:/DIRT_X/Gepatsch/Coupling/spotpy_cal/wasim_coupling/' 
    glaccells_grd      = 'glcalm%04i.grd'
    glaccodes_grd      = 'glidalm%04i.grd'
    glmb_grd           = 'glmbalm.grd'
    glmbold_grd        = 'glmb_oldalm.grd'
    station_area       = 57.6
    subbasins_areas    = {1:3.86, 
                        2:9.46, 201:1.68, 202:0.93,
                        3:4.11, 
                        4:6.71, 
                        5:11.59, 
                        6:4.89, 
                        7:36.98, 701:12.41, 702:4.14, 703:1.89, 704:2.98, 705:1.33, 706:1.01,
                        8:7.52, 801:14.16, 802:2.36, 803:2.69, 
                        9:16.58,
                        10:23.41,
                        11:23.38, 1101:6.42, 1102:26.02, 1103:1.78,
                        12:49.06}
    obs_runoff         = pd.read_csv(r'E:\DIRT_X\Gepatsch\Coupling\obs_discharge\gepatschalm_discharge.csv', index_col=0, parse_dates=True)['Q']
    obs_glmb_aux       = pd.read_csv(r'E:\DIRT_X\Gepatsch\Coupling\oggm_mass_balances\v1.6\MB_RGI60-11.00746.csv', index_col='time', parse_dates=True)
    obs_glmb           = obs_glmb_aux.resample('AS-SEP', origin='end').last()['mb_mm']

    target_id_runoff   = '11'
    target_id_glac     = '1102'
    sim_runoff         = 'qgkoalm'
    sim_glmb           = 'glmb2alm'
 
    calib_algorithm        = 'sceua'
    calib_runs             = 3000
        
    t1_cal                 = '1985-01-01 00:00:00'  
    t2_cal                 = '1998-12-31 00:00:00' 
