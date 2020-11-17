# -*- coding: utf-8 -*-
"""
Created on Thur July 13 18:43:31 2017
Correction of CoG movement
@author: Bas Blank
"""

import sys
import os

sys.path=['','/usr/local/lib64/python2.7/site-packages/distribute-0.6.28-py2.7.egg','/usr/local/lib64/python2.7/site-packages/matplotlib-0_unknown-py2.7-linux-x86_64.egg','/usr/local/lib64/python2.7/site-packages/pyparsing-2.0.6-py2.7.egg','/usr/local/lib64/python2.7/site-packages/cycler-0.9.0-py2.7.egg','/usr/local/lib64/python2.7/site-packages/pytz-2015.7-py2.7.egg','/usr/local/lib64/python2.7/site-packages/functools32-3.2.3_2-py2.7.egg','/usr/local/lib64/python2.7/site-packages/python_dateutil-2.4.2-py2.7.egg','/usr/local/lib64/python2.7/site-packages/six-1.10.0-py2.7.egg','/usr/lib/python27.zip','/usr/lib64/python2.7','/usr/lib64/python2.7/plat-linux2','/usr/lib64/python2.7/lib-tk','/usr/lib64/python2.7/lib-old','/usr/lib64/python2.7/lib-dynload','/usr/lib64/python2.7/site-packages','/usr/local/lib64/python2.7/site-packages','/usr/local/lib/python2.7/site-packages','/usr/lib/python2.7/site-packages','/usr/lib/python2.7/site-packages/setuptools-0.6c11-py2.7.egg-info','/usr/lib64/python2.7/site-packages/wx-2.8-gtk2-unicode']



import numpy as  np
print np.version.version
#import scipy.special as SP
#import scipy.misc as Misc
import scipy as sp
pi=np.pi
from scipy.interpolate import griddata
#from scipy.interpolate import Rbf


import spherepy as sph

if __name__=='__main__':
    Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
    Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
    sys.path.insert(0,Dir2)
    sys.path.insert(1,Dir)
    import Model_data2_top
    from Model_data2_top import *
    import Model_data2
    from Model_data2 import *
    os.chdir(complete_dir)
    print "The folder for Cog_correction is: "
    print os.getcwd()
       
    ####
            
    r_cog_vec_step = np.zeros([N_step,3])
    fac_load = np.zeros([2,3])
    
    if SLE_true  == 1: 
        Sea_lvl_allTimes_ice_eq = np.loadtxt('Data_in_Slvl.dat')*Density_water/Density_Ice
    else:
        Sea_lvl_allTimes_ice_eq = 0
    
    for l in range(2):
        
        fac_load[l,:]=4*pi*G/(2*l+1)*Radius[N_layer-1]
            
    for i in range(N_step):      
        
        load_data       = DataIce[(grid1+1)*i:(grid1+1)*(i+1),:]
        if SLE_true  == 1: 
            load_data_ocean = Sea_lvl_allTimes_ice_eq[(grid1+1)*i:(grid1+1)*(i+1),:]
            DATA_Load_ocean=sph.ScalarPatternUniform(load_data_ocean[:,0:-1])
            clm_seaLoad=sph.spht(DATA_Load_ocean,Degree,Degree) 
        else:
            DATA_Load_ocean=sph.ScalarPatternUniform(np.zeros([grid1+1,2*grid1]))
            clm_seaLoad=sph.spht(DATA_Load_ocean,Degree,Degree) 
            
        DATA_Load=sph.ScalarPatternUniform(load_data[:,0:-1])
        clm_iceLoad=sph.spht(DATA_Load,Degree,Degree)
    
        phi_load = (clm_iceLoad[1,:]*Density_Ice+clm_seaLoad[1,:]*Density_water)*fac_load[1,:]
        CoG_movement_lm = np.sqrt(3.0/np.pi)*phi_load/Gacc[N_layer-1]
        CoG_movement_lm_vector = np.zeros([3])
        CoG_movement_lm_vector[0] = 1.0/np.sqrt(2.0)*np.imag(CoG_movement_lm[2])
        CoG_movement_lm_vector[1] = 1.0/2.0*np.real(CoG_movement_lm[1])
        CoG_movement_lm_vector[2] = -1.0/np.sqrt(2.0)*np.real(CoG_movement_lm[2])
    
        r_cog_vec_step[i,:] = CoG_movement_lm_vector
        
    fileout=open('CoG_vector.dat','w')
    np.savetxt(fileout, r_cog_vec_step)
    fileout.close()
