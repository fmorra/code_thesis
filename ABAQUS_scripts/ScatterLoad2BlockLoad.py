import sys
import os

sys.path=['','/usr/local/lib64/python2.7/site-packages/distribute-0.6.28-py2.7.egg','/usr/local/lib64/python2.7/site-packages/matplotlib-0_unknown-py2.7-linux-x86_64.egg','/usr/local/lib64/python2.7/site-packages/pyparsing-2.0.6-py2.7.egg','/usr/local/lib64/python2.7/site-packages/cycler-0.9.0-py2.7.egg','/usr/local/lib64/python2.7/site-packages/pytz-2015.7-py2.7.egg','/usr/local/lib64/python2.7/site-packages/functools32-3.2.3_2-py2.7.egg','/usr/local/lib64/python2.7/site-packages/python_dateutil-2.4.2-py2.7.egg','/usr/local/lib64/python2.7/site-packages/six-1.10.0-py2.7.egg','/usr/lib/python27.zip','/usr/lib64/python2.7','/usr/lib64/python2.7/plat-linux2','/usr/lib64/python2.7/lib-tk','/usr/lib64/python2.7/lib-old','/usr/lib64/python2.7/lib-dynload','/usr/lib64/python2.7/site-packages','/usr/local/lib64/python2.7/site-packages','/usr/local/lib/python2.7/site-packages','/usr/lib/python2.7/site-packages','/usr/lib/python2.7/site-packages/setuptools-0.6c11-py2.7.egg-info','/usr/lib64/python2.7/site-packages/wx-2.8-gtk2-unicode']

Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
sys.path.insert(0,Dir2)
sys.path.insert(1,Dir)
import Model_data2_top
from Model_data2_top import *
os.chdir(complete_dir)
##print "The folder for ScatterLoad2BlockLoad is: "
##print os.getcwd()
import numpy as  np
print np.version.version
#import scipy.special as SP
#import scipy.misc as Misc
import scipy as sp
pi=np.pi
from scipy.interpolate import griddata
#from scipy.interpolate import Rbf


import spherepy as sph

def block_gen1(nodes_xyz,data,grid,Axis):
#   Generate n*2n grid and grid values, data at middle points
    a=nodes_xyz
    a2=data;
    if len(a2.shape)==1:
        a2=a2.reshape(len(a2),1)
        D_dim=1
    else:
        D_dim=a2.shape[1]
    
    if Axis=='Z':
        c=np.zeros(a.shape)
        c[:,0]=(a[:,0]**2+a[:,1]**2+a[:,2]**2)**0.5
        c[:,1]=sp.arccos(a[:,2]/c[:,0])
        c[:,2]=sp.arctan2(a[:,1],a[:,0])
    elif Axis=='X':
        c=np.zeros(a.shape)
        c[:,0]=(a[:,0]**2+a[:,1]**2+a[:,2]**2)**0.5
        c[:,1]=sp.arccos(a[:,0]/c[:,0])
        c[:,2]=sp.arctan2(a[:,2],a[:,1])
    elif Axis=='Y':
        c=np.zeros(a.shape)
        c[:,0]=(a[:,0]**2+a[:,1]**2+a[:,2]**2)**0.5
        c[:,1]=sp.arccos(a[:,1]/c[:,0])
        c[:,2]=sp.arctan2(a[:,0],a[:,2])

#    Mean_r=c[:,0].mean() #radius of the surface
    PoleN=np.nonzero(c[:,1]==c[:,1].min())[0]
    PoleS=np.nonzero(c[:,1]==c[:,1].max())[0]
    t1=np.nonzero(abs(c[:,2]-pi)<1e-4)
    t2=np.nonzero(abs(c[:,2]+pi)<1e-4)
    Meri=np.hstack([t1,t2])
#    DT=pi/grid

    grid_x, grid_y = np.mgrid[0:pi:(grid+1)*1j, -pi:pi:(grid*2+1)*1j]

    #padding northpole
    ex_c=np.zeros([2*grid,3])      
    ex_c2=np.linspace(-pi,pi,num=2*grid)
       
    ex_c[:,2]=ex_c2    
    c=np.vstack([c,ex_c])
    #padding southpole
    ex_c[:,1]=pi
    c=np.vstack([c,ex_c])
    
    ex_a2=np.zeros([2*grid,D_dim])
    ex_a2[:]=a2[PoleN,:]
    a2=np.vstack([a2,ex_a2])
    
    ex_a2[:]=a2[PoleS,:]
    a2=np.vstack([a2,ex_a2])
#    
#    ##padding -pi azimuth
    ex_c=np.zeros([Meri[0].size,3])
    ex_c[:,2]=-c[Meri,2]
    ex_c[:,1]=c[Meri,1]
    c=np.vstack([c,ex_c])
    
    ex_a2=np.zeros([Meri[0].size,D_dim])  
    ex_a2[:]=a2[Meri,:]
    a2=np.vstack([a2,ex_a2])


    grid_z = griddata(c[:,1:], a2[:,0], (grid_x, grid_y), method='linear')    

#    rbf = Rbf(c[:,1], c[:,2:], a2[:,0], epsilon=2)
#    grid_z= rbf(grid_x, grid_y)
    return grid_z#,grid_x,grid_y,

if __name__=='__main__':
##    if 'Model_data2_top' in (sys.modules.keys()):
##        reload(Model_data2_top)
##        from Model_data2_top import *
##    else:
##        Dir=r"/home/fabri/Earth_model_abaqus_SLE0"
##        os.chdir(Dir)
##        import Model_data2_top
##        from Model_data2_top import *
    Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
    Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
    sys.path.insert(0,Dir2)
    sys.path.insert(1,Dir)
    import Model_data2_top
    from Model_data2_top import *
    os.chdir(complete_dir)
    print "The folder for ScatterLoad2BlockLoad is: "
    print os.getcwd()
            
    DataIce = np.zeros([(grid1+1)*N_step,grid1*2+1])
    for i in range(N_step):      
        iceLoadName = 'IceLoad_' + Model_name + str(i+1) +'_field'
        fileIn = IceDir+'ice_load_xyzf_'+str(i+1)+'.txt'
        DataIce_scat=np.loadtxt(fileIn, delimiter=',') #kg/m^2    

        block_res=block_gen1(DataIce_scat[:,0:3],DataIce_scat[:,-1],grid1,Axis)
        block_res_height = block_res/Density_Ice
        DataIce[(grid1+1)*i:(grid1+1)*(i+1),:] = block_res_height
    np.savetxt('DataIce_in_int.dat',DataIce)
    
    ####
