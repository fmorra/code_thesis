## All constants and model data required to build the model
## Used both by the generator and iterator
## Author: Haiyang Hu
## Haiyang Hu: Build module (19-05-2014)
## Bas Blank: Added density of water (13-6-2016)

import numpy as np
import os
import sys

#if 'Model_data2_top' in (sys.modules.keys()):
#    reload(Model_data2_top)
#    from Model_data2_top import *
#else:
Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
sys.path.insert(0,Dir2)
sys.path.insert(1,Dir)
import Model_data2_top
from Model_data2_top import *
print ("The folder for Model_data2 is: ")
print (os.getcwd())
        
alpha = PolarCirkelRadiusDeg*Deg2Rad
alpha_north = N_PolarCirkelRadiusDeg*Deg2Rad

Data=np.loadtxt(Dir+'Layer data.txt')  
if Data.size==4:
    N_layer=1  # Number of layers
    
    Radius=np.array([Data[0]])
    Density=np.array([Data[1]])
    Density=np.hstack([Density,0,0])
    
    Y_mod=np.array([Data[2]])*(1+Poi)*2
    Vis=np.array([Data[3]])    
else:    
    N_layer=Data.shape[0]  # Number of layers
    
    Radius=Data[:,0]
    Density=Data[:,1]
    
    Density=np.hstack([Density,0])
    if Poi >=0.495:
        Y_mod=Data[:,2]*3
    else:
        Y_mod=Data[:,2]*(1+Poi)*2
    Vis=Data[:,3]

##if Vis[0]==0:
##    flag_lc=1
Density0=np.zeros(N_layer)
Gacc=np.zeros(N_layer) # Gravitational acceleration
for i in range(N_layer):
    Gacc[i]=(4.0/3*G*pi*Radius[0]**3*Density[0]+np.sum([4.0/3*G*pi*Density[j+1]*(Radius[j+1]**3-Radius[j]**3) for j in range(i) ]))/Radius[i]**2
    if i == 0:
        Density0[i] = Density[0]
    else:
        Density0[i] = Radius[i-1]**3/Radius[i]**3*(Density0[i-1]-Density[i])+Density[i]
   
### Keep grid consistent
if sedimentation_history == 1:
    DataSediments=np.loadtxt(Dir + 'DataSediments_in.dat')
if scatterGrid_load==0:
    
    DataIce=np.loadtxt(Dir + 'DataIce_in.dat')
    grid1=DataIce.shape[1]/2
    if initial_ice == 1:
        Initial_Ice = np.loadtxt(Dir + 'DataIceInitial_in.dat')
    if filter_elastic==1:
        DataIce[(grid1+1)*(N_step-1):(grid1+1)*N_step] = DataIce[(grid1+1)*(N_step-2):(grid1+1)*(N_step-1)]
        
    if SLE_true == 1 :##and Only_Eustatic == 0: 
        Topo_present = np.loadtxt(Dir + 'Topo_present.dat')
        if initial_ice==1:
            Topo_present = Topo_present+Initial_Ice*Density_Ice/Density_water
        Topo_present_matrix = np.tile(Topo_present,(N_step,1))
        Ocean_Function = 1*(Topo_present_matrix<0)
    resolution = 360.0/(2*grid1)
else:
    print 'Start conversion from scattered gridded load to equi-angular gridded load\n'
    CHECK=os.system('python2.7 /home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2/ScatterLoad2BlockLoad.py') ##1
    print 'Finish initial sub-Python process\n'
    if CHECK!=0:
        print('Error')
        
    DataIce=np.loadtxt('DataIce_in_int.dat')
    SLE_true = 0   
    
##if initial_ice == 1:
##    Initial_Ice = np.loadtxt(Dir + 'DataIceInitial_in.dat')
##if SLE_true == 1 :##and Only_Eustatic == 0: 
##    Topo_present = np.loadtxt(Dir + 'Topo_present.dat')
##    if initial_ice==1:
##        Topo_present = Topo_present+Initial_Ice*Density_Ice/Density_water
##    Topo_present_matrix = np.tile(Topo_present,(N_step,1))
##    Ocean_Function = 1*(Topo_present_matrix<0)
##if sedimentation_history == 1:
##    DataSediments=np.loadtxt(Dir + 'DataSediments_in.dat')
##if scatterGrid_load==0:   
##    DataIce=np.loadtxt(Dir + 'DataIce_in.dat')
##    grid1=DataIce.shape[1]/2
##
##    if filter_elastic==1:
##        DataIce[(grid1+1)*(N_step-1):(grid1+1)*N_step] = DataIce[(grid1+1)*(N_step-2):(grid1+1)*(N_step-1)]
##        
##    resolution = 360.0/(2*grid1)
##else:
##    print 'Start conversion from scattered gridded load to equi-angular gridded load\n'
##    CHECK=os.system('python2.7 ' + Dir + 'ScatterLoad2BlockLoad.py') ##1
##    print 'Finish initial sub-Python process\n'
##    if CHECK!=0:
##        keyboarddd
##        
##    DataIce=np.loadtxt('DataIce_in_int.dat')
##    grid1=DataIce.shape[1]/2

##elif SLE_true == 1 and Only_Eustatic == 1: 
##    Raw_Ocean_Function = np.loadtxt('OceanFunction.dat')
##
##    #matrix operation for if (ice==true && ocean==true) ->ocean = false
##    np.seterr(all='ignore')
##    Normalized_Ice = DataIce/DataIce
##    Where_Nans = np.isnan(Normalized_Ice)
##    Normalized_Ice[Where_Nans] =0
##    Ocean_Function = Raw_Ocean_Function-np.multiply(Normalized_Ice,Raw_Ocean_Function)


### Radius of polar cirkel

R_PolCir = PolarCirkelRadiusDeg * Radius[N_layer-1]
