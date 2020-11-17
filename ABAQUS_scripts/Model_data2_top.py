import numpy as np
import os
import sys

## Working directory
Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
sys.path.insert(0,Dir2)
sys.path.insert(1,Dir)
run_number = 21
dir_name = "results_run_" + str(run_number)
complete_dir = os.path.join(Dir,dir_name)
if not os.path.exists(complete_dir):
    os.mkdir(complete_dir)
os.chdir(complete_dir)
IceDir=Dir+"IceFiles/"
CPUs = 16 #Number of CPU's allocated

## Constants
G=6.6732e-11 #85e-11
pi=np.pi
ka=1000*365.25*24*3600
Density_water = 1000.0 #999.972
Density_Ice = 931.0
Density_Sed = 2300.0
Deg2Rad = pi/180

Omega=0#pi/24/3600
## Layer data are defined in file: Layer data.txt in the working directory

flag_lc=1  # Indicator for liquid core, 1 for liquid (core not in the simulation), 0 for solid. 


### Set variable for centrifugal effect (set 1 for true, set 0 for false)
Centrifugal_on = 0
### fix centre of the model
fixed_point =1

### Set variable for simple loading without ice mass changes (set 1 for true, set 0 for false)
SLE_true = 1
Only_Eustatic = 0
initial_ice = 1
sedimentation_history = 0

if SLE_true == 0:
    No_IceLoss = 1
else:
    No_IceLoss = 0
    
Ice_sealvl = Only_Eustatic ### sea level rise only determined by ice equivelent mass loss
###

## Model name
Model_name='Earth'

CETOL=1e-2
    
## Time step 
##Time=np.array([0.01,0.1,0.2,0.3,0.5,0.8])*ka # Time steps
##Time=np.array([0.test01,0.1,0.2,0.3,0.5,0.8,1.1,1.5,2,3,4,6,10])*ka # Time steps 13step
##Time=np.array([0.01,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.65,0.8,1,1.3,1.6,2,3,4,5,6,8,10])*ka # Double Time steps
##Time=np.array([0.01,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.65,0.8,1,1.3,1.6,2,3,4,5,6,8,10,12,15])*ka # Double+ Time steps
##Time=np.array([0.01,0.1,0.2,0.3,0.5,0.8,1.1,1.5,2,3,4,6,10,15])*ka # Time steps 14step
##Time=np.array([0.01,1,2,3,4,6,8,10])*ka 
##Time=np.array([10])*ka # Time steps 13step
##Time=np.hstack([np.linspace(0,1,15)[0:-1], np.linspace(1,2,7)[0:-1], 2,3,5,7,10])*ka 
##Time=np.array([1])*ka 
##Time=np.array([0.001,10.0])*ka
##Time=np.array([90.0,91.0])*ka
##Time=np.array([90.0,91.0, 92.0])*ka
##Time=np.array([90.0, 91.0, 91.5, 95.0, 95.5, 100.0, 100.5])*ka

##Time = np.array([0.0001, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 ])*ka
##Time = np.array([1.0, 31.0, 31.051, 31.102, 31.103, 31.104, 31.105, 31.106, 31.107, 31.108, 31.109, 31.110, 31.111, 31.112, 31.113, 31.114, 31.115  ])*ka
Time = np.array([1.0, 31.0, 31.051, 31.102, 31.104, 31.106, 31.108, 31.110, 31.112, 31.114, 31.115])*ka
##Time = np.array([102.0, 104.5, 107.0, 108.0, 109.0, 110.5, 112.0, 113.5, 115.0, 117.0, 118.0, 119.0, 120.0, 121.0, 121.5, 121.9, 122.002, 122.005, 122.008, 122.011, 122.014, 122.015, 122.016])*ka
##Time = np.array([10.0, 50.0, 90.0, 110.0, 112.5, 115.0, 116.5, 118.0, 119.5, 121.0, 123.0, 124.0, 126.0, 127.0, 128.0, 129.0, 130.0, 130.01])*ka

##Time[0]=0.001*ka
N_step=Time.shape[0]

## Rampload
Rampload_enabled=1
filter_elastic=0
scatterGrid_load=1
## if scatterGrid is set to 1, define a grid_resolution [pixels/deg] (dont go below 1). For normal loading these parameters will be determined based on the shape of the input
resolution=4
grid1= 180*resolution
Axis='Y'

## Spherical harmonics degree limit, seeds of FEM nodes distance
Degree=360
DisableSH_0_1 = 0 #(1 = yes| 0 = no)
Compensate_cog = 1#!!!!!!!!!!!!!! CHECK!!!!!!!!!!!!!!!
Sphere_int=0
initialize_interpol_lists=0

Seeds=200e3
PolarCirkelRadiusDeg = 15.0 #Use doubles!
N_PolarCirkelRadiusDeg = 4.0
meshtickness = 66#40
Hemisphere_seeding_bias = 8
Polar_seeding_bias = 5
R_Target_seed = 50E3#50E3
Plane_Target_seed = 25E3#42E3
LM_res = 1
hires_rotate_grid = 1

Coarse_Layers = 1
## Pre-Processing
##Poi=0.5-1e-22
##Poi=0.495
Poi=0.49999999
##Poi=0.28
## ROTATE LOCATION OF HIGH-RESOLUTION (hr) MESH AREA TO DIFFERENT FOCAL POINT



# Set 0 to leave hi-res mesh area at South Pole

# Set 1 to rotate hi-res mesh area to any other location

if hires_rotate_grid == 1:
#Define centre of new hi-res mesh area
    phi_hr = -76.0#58.0 #degrees latitude
    theta_hr = -108.3#11.0 #degrees longitude.
