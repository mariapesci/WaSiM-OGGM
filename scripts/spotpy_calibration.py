import glob
import os
from datetime import datetime
import pandas as pd
import numpy as np
import copy
import sys
import spotpy 
from spotpy.parameter import Uniform
#from spotpy_settings import wasim_model, wasim_settings
from spotpy_coupling_settings import wasim_coupling_model, coupling_settings


ws = coupling_settings

sys.path.append(".")


def adjust_series(y_obs, y_sim, t_ini=None, t_fin=None):
    
    if t_ini or t_fin is None:
        ind_sim = y_sim.index
        t1s = ind_sim[0]
        t2s = ind_sim[-1]
        
        ind_obs = y_obs.index
        t1o  = ind_obs[0]
        t2o  = ind_obs[-1]
        
        t1 = max(t1o,t1s)
        t2 = min(t2o,t2s)
    else:
        t1 = t_ini
        t2 = t_fin

    y_b = y_obs[(y_obs.index>=t1) & (y_obs.index<=t2)]
    y_s = y_sim[(y_sim.index>=t1) & (y_sim.index<=t2)]

    return y_b, y_s  
    
def bench_nse(y_obs, y_sim, t_ini=None, t_fin=None):

    y_bench = (y_obs.groupby('{:%m-%d}'.format).mean()).to_dict()
    y_o = pd.DataFrame(copy.deepcopy(y_obs))
    y_o['month_day'] = y_o.index.strftime('%m-%d')
    y_o['bench'] = 0.
    for i in range(len(y_o['month_day'])):
        if y_o['month_day'][i] in y_bench:
            y_o['bench'][i] = y_bench[y_o['month_day'][i]]
    y_b, y_s = adjust_series(y_o[y_o.columns[0]], y_sim, t_ini, t_fin)      
    y_ben,_ = adjust_series(y_o['bench'], y_sim, t_ini, t_fin)
       
    be = 1-(np.sum((y_b-y_s)**2)/np.sum((y_b-y_ben)**2))

    return be


class spotpy_setup():
    
    def __init__(self, obs_runoff, obs_glmb, t1, t2, add_t0=True, add_t0r=False, add_tros=False, add_trs=False,
                 add_mf=True, add_ddf=False, add_snow_param=True, add_ice_param=True,
                 add_kd=True, add_ki=True, add_q0=True):
        self.params = []
        if add_t0:
            self.params.append(spotpy.parameter.Uniform('t0', -1.0, 2.0, optguess=0.0))
        if add_t0r:
            self.params.append(spotpy.parameter.Uniform('t0r', -1.0, 2.0, optguess=-0.4))
        if add_tros:
            self.params.append(spotpy.parameter.Uniform('tros', 0.5, 1.0, optguess=0.5))
        if add_trs:
            self.params.append(spotpy.parameter.Uniform('trs', -1.0, 2.0, optguess=0.0))
        if add_kd:        
            self.params.append(spotpy.parameter.Uniform('kd', 5.0, 300., optguess=270.))
        if add_ki:        
            self.params.append(spotpy.parameter.Uniform('ki', 5.0, 300., optguess=210.))
        if add_q0:        
            self.params.append(spotpy.parameter.Uniform('q0', 0.01, 1.0, optguess=1.0))
        if add_mf:        
            self.params.append(spotpy.parameter.Uniform('mf', 1.2, 4.0, optguess=2.0))
        if add_ddf:
            self.params.append(spotpy.parameter.Uniform('DDFi', 2., 10., optguess=4.0))
            self.params.append(spotpy.parameter.Uniform('DDFf', 2., 8., optguess=3.5))
            self.params.append(spotpy.parameter.Uniform('DDFs', 2., 6., optguess=2.0))
        # if add_ddf:
        #     self.params.append(spotpy.parameter.Uniform('ki', 1., 20., optguess=20.0))
        #     self.params.append(spotpy.parameter.Uniform('kf', 100., 1000., optguess=500.0))
        #     self.params.append(spotpy.parameter.Uniform('ks', 10., 100., optguess=80.0))
        if add_snow_param:
            self.params.append(spotpy.parameter.Uniform('min_slope', 0., 90., optguess=60.0))
            #self.params.append(spotpy.parameter.Uniform('fr_snow', 0., 1., optguess=0.076))
            #self.params.append(spotpy.parameter.Uniform('lwin', 0.8, 1.2, optguess=1.15))
            self.params.append(spotpy.parameter.Uniform('lwout', 0.8, 1.2, optguess=1.07))
        if add_ice_param:
            self.params.append(spotpy.parameter.Uniform('k_ice', 1, 20., optguess=20.0))
            self.params.append(spotpy.parameter.Uniform('k_firn', 100., 1000., optguess=500.0))
            self.params.append(spotpy.parameter.Uniform('k_snow', 10., 100., optguess=80.0))     
            #self.params.append(spotpy.parameter.Uniform('firn_stack', 7., 10., optguess=10.0))      
            #self.params.append(spotpy.parameter.Uniform('rate_we', 1., 2., optguess=1.5))  
        
        
        
        self.eval_runoff = obs_runoff
        self.eval_glmb = obs_glmb
        self.t1_cal = t1
        self.t2_cal = t2
        #for ii,pi in enumerate(self.parameters()['name']):
            #print('%i ... %s' % (ii, pi))
       
        
    def parameters(self):
        
        return spotpy.parameter.generate(self.params)
    
    
    def simulation(self, vector):

        param_var = list(self.parameters()['name'])
        
        for pi, pj in zip(param_var, vector):
            print('%s = %s' % (pi, pj))   
        
        # Running only wasim
        # model = wasim_model(ws.working_dir, ws.ctrl_master, ws.ctrl_cal, 
        #                                   param_var, vector, ws.sim_runoff, ws.sim_glmb, 
        #                                   ws.target_id_runoff, ws.target_id_glac, 
        #                                   res_include_stats=False, verbose=True)
        
        # Running coupling model
        model = wasim_coupling_model(ws, ws.ctrl_master, ws.ctrl_cal, 
                                          param_var, vector, ws.sim_runoff, ws.sim_glmb, 
                                          ws.target_id_runoff, ws.target_id_glac)

        self.runoff, self.glmb_annual = model.run()

        return self.runoff, self.glmb_annual

   
    def evaluation(self):
        
        return self.eval_runoff, self.eval_glmb

    
   
    def objectivefunction(self, simulation, evaluation, params=None):
 
        sim1 = self.runoff
        sim2 = self.glmb_annual
        obs1 = self.eval_runoff
        obs2 = self.eval_glmb

        runoff_obs, runoff_sim = adjust_series(obs1, sim1, self.t1_cal, self.t2_cal)
        glmb_obs, glmb_sim = adjust_series(obs2, sim2, self.t1_cal, self.t2_cal)
                
        obj1 = spotpy.objectivefunctions.kge(runoff_obs, runoff_sim)
        obj2 = bench_nse(obs1, sim1, self.t1_cal, self.t2_cal)
        obj3 = spotpy.objectivefunctions.bias(runoff_obs, runoff_sim)
        obj4 = spotpy.objectivefunctions.rsr(glmb_obs, glmb_sim)
        
        # w1, w2, w3 = 0.63, 0.07, 0.3
        # multi_obj = w1*(1-obj1)+(1-w1)*(1-obj2)+w2*abs(obj3)+w3*obj4
        w1, w2, w3, w4 = 0.23, 0.4, 0.07, 0.3
        multi_obj = w1*(1-obj1)+w2*(1-obj2)+w3*abs(obj3)+w4*obj4
        
        return multi_obj # [obj1, obj2]
 

model_run = spotpy_setup(ws.obs_runoff, ws.obs_glmb, ws.t1_cal, ws.t2_cal)
#(working_dir, ctrl_master, ctrl_new, param_var, param_val, obs_runoff, obs_glmb, sim_runoff, sim_glmb)
sampler = spotpy.algorithms.sceua(model_run, dbformat='csv', dbname='cal_coup_init_1969')
sampler.sample(ws.calib_runs, ngs=7, kstop=3, peps=0.1, pcento=0.1)
# ngs: number of complexes in the initial population
# kstop: number of shuffling loops in which the criterion value must change by the given % before optimiz. is terminated
# peps: convergence level for parameter set
# pcento: % by which the criterion value must change in a given kstop of shuffling loops to continue optimization

#results_spotpy = spotpy.analyser.load_csv_results('cal_coup_3012', usecols=0)
# results_spotpy = np.genfromtxt(r'E:\DIRT_X\Gepatsch\Coupling\spotpy_cal\wasim_coupling\cal_coup_0301_2.csv', delimiter=',', usecols=0)
# import matplotlib.pyplot as plt
# fig= plt.figure(1,figsize=(9,5))
# plt.plot(results_spotpy)#['like1'])
# plt.show()
# plt.ylabel('Multi objective function')
# plt.xlabel('Iteration')
