# -*- coding: utf-8 -*-
"""
Created on Tue May 20 11:51:41 2014
Post-processing for the model
@authors: Bas Blank
          Haiyang Hu
         

28-10-2016 hh: 1. Replace scipy.Special.lpmv with spherePy package
               2. block_gen1 interpolate the data in the corner points 
"""
# This is another version of the code used to save results at each timestep of the sea level equation
import sys
import os
import shutil

sys.path=['','/usr/local/lib64/python2.7/site-packages/distribute-0.6.28-py2.7.egg','/usr/local/lib64/python2.7/site-packages/matplotlib-0_unknown-py2.7-linux-x86_64.egg','/usr/local/lib64/python2.7/site-packages/pyparsing-2.0.6-py2.7.egg','/usr/local/lib64/python2.7/site-packages/cycler-0.9.0-py2.7.egg','/usr/local/lib64/python2.7/site-packages/pytz-2015.7-py2.7.egg','/usr/local/lib64/python2.7/site-packages/functools32-3.2.3_2-py2.7.egg','/usr/local/lib64/python2.7/site-packages/python_dateutil-2.4.2-py2.7.egg','/usr/local/lib64/python2.7/site-packages/six-1.10.0-py2.7.egg','/usr/lib/python27.zip','/usr/lib64/python2.7','/usr/lib64/python2.7/plat-linux2','/usr/lib64/python2.7/lib-tk','/usr/lib64/python2.7/lib-old','/usr/lib64/python2.7/lib-dynload','/usr/lib64/python2.7/site-packages','/usr/local/lib64/python2.7/site-packages','/usr/local/lib/python2.7/site-packages','/usr/lib/python2.7/site-packages','/usr/lib/python2.7/site-packages/setuptools-0.6c11-py2.7.egg-info','/usr/lib64/python2.7/site-packages/wx-2.8-gtk2-unicode']
Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
sys.path.insert(0,Dir2)
sys.path.insert(1,Dir)

import Model_data2_top
from Model_data2_top import *
import Model_data2
from Model_data2 import *

import numpy as  np
print (np.version.version)
#import scipy.special as SP
#import scipy.misc as Misc
import scipy as sp
pi=np.pi
from scipy.interpolate import griddata
from numpy import complex64

os.chdir(complete_dir)
print ("The current directory for sph_tools is:\n")
print (os.getcwd())
#from scipy.interpolate import Rbf

import spherepy as sph
def block_gen_sphere_int(data,grid,layer):
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

    Mean_r=c[:,0].mean() #radius of the surface
    PoleN=np.nonzero(c[:,1]==c[:,1].min())[0]
    PoleS=np.nonzero(c[:,1]==c[:,1].max())[0]
    t1=np.nonzero(abs(c[:,2]-pi)<1e-4)
    t2=np.nonzero(abs(c[:,2]+pi)<1e-4)
    Meri=np.hstack([t1,t2])
#    DT=pi/grid

#    grid_x, grid_y = np.mgrid[DT/2:pi-DT/2:grid*1j, -pi+DT/2:pi-DT/2:grid*2j]
    grid_x, grid_y = np.mgrid[0:pi:(grid+1)*1j, -pi:pi:(grid*2+1)*1j]

    #padding northpole
    ex_c=np.zeros([2*grid,3])      
#    ex_c2=np.linspace(-pi+DT/2,pi-DT/2,num=2*grid)
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
    return grid_z,Mean_r#,grid_x,grid_y,

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

    Mean_r=c[:,0].mean() #radius of the surface
    PoleN=np.nonzero(c[:,1]==c[:,1].min())[0]
    PoleS=np.nonzero(c[:,1]==c[:,1].max())[0]
    t1=np.nonzero(abs(c[:,2]-pi)<1e-4)
    t2=np.nonzero(abs(c[:,2]+pi)<1e-4)
    Meri=np.hstack([t1,t2])
#    DT=pi/grid

#    grid_x, grid_y = np.mgrid[DT/2:pi-DT/2:grid*1j, -pi+DT/2:pi-DT/2:grid*2j]
    grid_x, grid_y = np.mgrid[0:pi:(grid+1)*1j, -pi:pi:(grid*2+1)*1j]

    #padding northpole
    ex_c=np.zeros([2*grid,3])      
#    ex_c2=np.linspace(-pi+DT/2,pi-DT/2,num=2*grid)
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
    return grid_z,Mean_r#,grid_x,grid_y,

def MI (grid_Z,R):  ##need to change
##  calculate change of moment of inertia due to deformation
## change from v4: calculated from grid being defined on corners 
    I=np.zeros((3,3))
    grid=grid_Z.shape[0]-1
    Theta=np.linspace(0,np.pi,grid+1)
    

    S=4*pi*R**2/(sum(abs(np.sin(Theta)))*(2*grid))
    
    weight=S*np.sin(Theta)
    weight=weight.reshape(grid+1,1)

    grid_t, grid_p = np.mgrid[0:pi:(grid+1)*1j, -pi:pi:(grid*2+1)*1j]
    grid_x=R*np.sin(grid_t)*np.cos(grid_p)
    grid_y=R*np.sin(grid_t)*np.sin(grid_p)
    grid_z=R*np.cos(grid_t)
    
    grid_x=grid_x[:,0:-1]
    grid_y=grid_y[:,0:-1]
    grid_z=grid_z[:,0:-1]
    
    I11=sum(sum((grid_y**2+grid_z**2)*grid_Z[:,0:-1]*weight))
    I22=sum(sum((grid_x**2+grid_z**2)*grid_Z[:,0:-1]*weight))
    I33=sum(sum((grid_y**2+grid_x**2)*grid_Z[:,0:-1]*weight))
    I12=sum(sum((-grid_y*grid_x)*grid_Z[:,0:-1]*weight))
    I13=sum(sum((-grid_z*grid_x)*grid_Z[:,0:-1]*weight))
    I23=sum(sum((-grid_z*grid_y)*grid_Z[:,0:-1]*weight))
    I=np.array([I11,I22,I33,I12,I13,I23])
    return I


    
class wrong_length(Exception):
    def __str__(self): return 'CLM,SLM,GRID must have the same number of elements'    
class m_1(Exception):
    def __str__(self): return 'm must be scalar'
    
class block_size(Exception):
    def __str__(self): return 'block must be size 2n*n'  

if __name__=='__main__':
##    Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
##    sys.path.insert(1,Dir)
##    import Model_data2
##    from Model_data2 import *
##    os.chdir(complete_dir)
#==============================================================================
#    print clm
#==============================================================================

##    load_data=np.loadtxt('Load_beforeSH_wouter.dat')
##            
##    DATA_Load=sph.ScalarPatternUniform(load_data)
##    clm= sph.spht(DATA_Load,Degree,Degree)
##    
##    clm[0,0]=0
##    clm[1,:]=[0,0,0]
##    
##    fileout=open('Int_count.dat','w')
##    clm_output = np.zeros([Degree+1,Degree*2+1])
##    for row in range(2,Degree+1):
##        integer = -1*row
##        fileout.write('coef clm =')
##        fileout.write('%s \n' % clm[row,-1])
##        fileout.write('%s \n' % clm[row,1])
##        for col in range((Degree-row),(Degree+row+1)):
##            fileout.write('Int=')
##            fileout.write('%s \n' % integer)
##            fileout.write('row+col \n')
##            fileout.write('%s \n' % row)
##            fileout.write('%s \n' % col)
##            clm_output[row, col] = clm[row,integer]
##            integer=integer+1
##
##    fileout.close()        
##    ##create output file3
##    fileout=open('Load_afterSH.dat','w')
##    np.savetxt(fileout,clm_output)
##    fileout.close()
#==============================================================================
#    N_step=1     ##Only major chnge from sph_tools_ult2.py
#==============================================================================     
    # Seems we do not need to change the folder by adding Dir as these files are 
    # saved together with the cae and odb in the subfolder

    Config=np.loadtxt('Nodes_config.dat')
    Config=np.array([int(i) for i in Config])
    Data=np.loadtxt('Data_output.dat')
    Data_Hor = np.loadtxt('Data_output_hor.dat')
    K2=np.zeros([1,N_step])
    print('test')
    
#     clm_stock3 = list((N_layer-1)*(Degree+1))

#    slm_stock=np.zeros(((Degree+1)*N_layer, Degree+1))
    Axis='Y'
    for i in range(N_layer):
        exec("Data_in"+str(i)+"=np.zeros([(grid1+1)*N_step,(grid1*2+1)])")
        exec("Deflection_Hor"+str(i)+"=np.zeros([(grid1+1)*N_step,(grid1*2+1)])")
    
        
    p20=np.array([])   
    clm_stock2 = np.zeros([N_layer*(Degree+1),2*Degree+1],dtype=complex)
    MoI=np.zeros([N_step,6])
    temp_moi=np.zeros([N_layer,6])

    # Allocate the dimension of the deflection matrices and their component matrices

    Deflection_Surface = np.real(np.zeros([grid1+1,grid1*2+1,N_step]))
    Deflection_surface = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
    geoid              = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
    Deflection = np.real(np.zeros([grid1+1,grid1*2+1]))
    
    dx_surface = np.real(np.zeros([grid1+1,grid1*2+1,N_step]))
    dy_surface = np.real(np.zeros([grid1+1,grid1*2+1,N_step]))          
    dz_surface = np.real(np.zeros([grid1+1,grid1*2+1,N_step]))
    
    r_cog_vec_Steps = np.zeros([N_step,3])

    # If we disable the spherical harmonics of degrees 0 and 1 we have to allocate some more files
    if DisableSH_0_1==1:
        Deflection_SurfaceSH01 = np.zeros([grid1+1,grid1*2+1,N_step])
        dx_surface_SH01 = np.zeros([grid1+1,grid1*2+1,N_step])
        dy_surface_SH01 = np.zeros([grid1+1,grid1*2+1,N_step])
        dz_surface_SH01 = np.zeros([grid1+1,grid1*2+1,N_step])    

    # Remember that n_step is the number of the ice history epochs
    # Define initial convergence value and extract iteration number
    
    Deflection_layers = np.real(np.zeros([(N_layer)*(grid1+1),grid1*2+1]))
    clm_stock=list()
    for j in range(N_layer):
        # Homogeneous Earth case
        if N_layer==1:
            ur=Data[i*Config:(i+1)*Config,:]
            dU=Data_Hor[i*Config:(i+1)*Config,:]
        # Layered Earth
        else:
            ur=Data[(i*sum(Config)+sum(Config[0:j])):(i*sum(Config)+sum(Config[0:j+1])),:]
            dU=Data_Hor[(i*sum(Config)+sum(Config[0:j])):(i*sum(Config)+sum(Config[0:j+1])),:]
        
        if ur.max()!=0:

            block_data1,r1=block_gen1(ur[:,0:3],ur[:,-1],grid1,Axis)
            block_data_dx,r1 = block_gen1(dU[:,0:3],dU[:,3],grid1,Axis)
            block_data_dy,r1 = block_gen1(dU[:,0:3],dU[:,4],grid1,Axis)
            block_data_dz,r1 = block_gen1(dU[:,0:3],dU[:,5],grid1,Axis)
            
            ## !input definition is [0 360] long while output definition is [-180 180]!
            data_reverse = np.zeros([grid1+1,grid1*2+1])
            data_reverse[:, :(grid1+1)] = block_data1[:,(grid1):]
            data_reverse[:, (grid1):] = block_data1[:,:(grid1+1)]
            block_data1 = data_reverse
            
            data_reverse = np.zeros([grid1+1,grid1*2+1])
            data_reverse[:, :grid1] = block_data_dx[:,(grid1+1):]
            data_reverse[:, (grid1+1):] = block_data_dx[:,:grid1]
            block_data_dx = data_reverse
            
            data_reverse = np.zeros([grid1+1,grid1*2+1])
            data_reverse[:, :grid1] = block_data_dy[:,(grid1+1):]
            data_reverse[:, (grid1+1):] = block_data_dy[:,:grid1]
            block_data_dy = data_reverse
            
            data_reverse = np.zeros([grid1+1,grid1*2+1])
            data_reverse[:, :grid1] = block_data_dz[:,(grid1+1):]
            data_reverse[:, (grid1+1):] = block_data_dz[:,:grid1]
            block_data_dz = data_reverse
            
            if j == (N_layer-1):
                Deflection_Surface[:,:,i] = block_data1
                dx_surface[:,:,i] = block_data_dx
                dy_surface[:,:,i] = block_data_dy
                dz_surface[:,:,i] = block_data_dz
                fileout=open('check_y_data.dat','w')
                np.savetxt(fileout,dy_surface[:,:,i])
            Deflection =  block_data1##delete entire else statement also delete printing code later
                    
        else:
            print ([i,j])
            block_data1=np.zeros([grid1+1,grid1*2+1])
            block_data_dx=np.zeros([grid1+1,grid1*2+1])  
            block_data_dy=np.zeros([grid1+1,grid1*2+1])  
            block_data_dz=np.zeros([grid1+1,grid1*2+1])  
        
            
        Deflection_layers[(j)*(grid1+1):(j+1)*(grid1+1),:] = Deflection
        if j == (N_layer-1):       
            eval("np.savetxt(os.path.join(complete_dir, 'Deflection_layers_step"+str(i)+".dat'),Deflection_layers)")
            
        DATA=sph.ScalarPatternUniform(block_data1[:,:-1])
        clm=sph.spht(DATA,Degree,Degree)
        
        DATA=sph.ScalarPatternUniform(block_data_dx[:,:-1])
        dx_lm=sph.spht(DATA,Degree,Degree)
        
        DATA=sph.ScalarPatternUniform(block_data_dy[:,:-1])
        dy_lm=sph.spht(DATA,Degree,Degree)
        
        DATA=sph.ScalarPatternUniform(block_data_dz[:,:-1])
        dz_lm=sph.spht(DATA,Degree,Degree)

##            clm[0,0]=0
##            clm[1,:]=[0,0,0]
##            
        clm_stock.append(clm)
#==============================================================================
#   Calculate Change of MOI of each layer 
       
        ## output surface deformation without 0 and 1 degree components 
        if DisableSH_0_1==1:            
            if j == (N_layer-1):
                clm_SH01 = clm[0:1,:]
                block_data_refined2 = np.real(sph.ispht(clm_SH01,grid1+1,grid1*2).cdata)
                block_data_refined2 = np.hstack([block_data_refined2, block_data_refined2[:,0:1]])
                 
                Deflection_SurfaceSH01[:,:,i] = block_data_refined2
                Deflection_Surface[:,:,i] = Deflection_Surface[:,:,i] - Deflection_SurfaceSH01[:,:,i]
                                
                dx_lm_SH01 = dx_lm[0:1,:]
                block_data_refined2 = np.real(sph.ispht(dx_lm_SH01,grid1+1,grid1*2).cdata)
                block_data_refined2 = np.hstack([block_data_refined2, block_data_refined2[:,0:1]])
                 
                dx_surface_SH01[:,:,i] = block_data_refined2
                dx_surface[:,:,i] = dx_surface[:,:,i]# - dx_surface_SH01[:,:,i]
                                
                dy_lm_SH01 = dy_lm[0:1,:]
                block_data_refined2 = np.real(sph.ispht(dy_lm_SH01,grid1+1,grid1*2).cdata)
                block_data_refined2 = np.hstack([block_data_refined2, block_data_refined2[:,0:1]])
                 
                dy_surface_SH01[:,:,i] = block_data_refined2
                dy_surface[:,:,i] = dy_surface[:,:,i]# - dy_surface_SH01[:,:,i]
                                
                dz_lm_SH01 = dz_lm[0:1,:]
                block_data_refined2 = np.real(sph.ispht(dz_lm_SH01,grid1+1,grid1*2).cdata)
                block_data_refined2 = np.hstack([block_data_refined2, block_data_refined2[:,0:1]])
                 
                dz_surface_SH01[:,:,i] = block_data_refined2
                dz_surface[:,:,i] = dz_surface[:,:,i]# - dz_surface_SH01[:,:,i]
                    
        block_data_refined=np.zeros(block_data1.shape)
##            clm[0,0]=0
##            clm[1,:]=[0,0,0]
##            clm_stock.append(clm)
        
        block_data_refined=np.real(sph.ispht(clm,grid1+1,grid1*2).cdata)
        block_data_refined=np.hstack([block_data_refined, block_data_refined[:,0:1]])
         
        Load_check=np.real(sph.ispht(clm,grid1+1,grid1*2).cdata)
        Load_check=np.hstack([Load_check, Load_check[:,0:1]])
         
            ##create output file3
        fileout=open('Load_afterSH_InverseTransform.dat','w')
        np.savetxt(fileout, Load_check)
        fileout.close()
         
        temp_moi[j,:]=MI(block_data_refined,Radius[j])*(Density[j+1]-Density[j])

#==============================================================================
#       For triggering loads
    Ice_data_allTimes = DataIce*Density_Ice
        
    if No_IceLoss  == 0:
        Sea_lvl_allTimes = np.loadtxt('Data_in_Slvl.dat')*Density_water
    else:
        Sea_lvl_allTimes = 0
        
    if sedimentation_history == 1:
        Sediments_allTimes = DataSediments*Density_Sed
    else:
        Sediments_allTimes = 0
        
    Total_load_allTimes =  Ice_data_allTimes +  Sea_lvl_allTimes + Sediments_allTimes
     
    if filter_elastic==1 and i==(N_step-1):
        load_data = Total_load_allTimes[(grid1+1)*(i-1):(grid1+1)*i,:]
    else:
        load_data = Total_load_allTimes[(grid1+1)*i:(grid1+1)*(i+1),:]
    ## !input definition is [0 360] long while output definition is [-180 180]!
##        load_data_reverse = np.zeros([grid1+1,grid1*2+1])
##        load_data_reverse[:, :grid1] = load_data[:,(grid1+1):]
##        load_data_reverse[:, (grid1+1):] = load_data[:,:grid1]
    load_data_reverse = load_data
##        fileout=open('reverse_check.dat','w')
##        np.savetxt(fileout, load_data_reverse)
##        fileout.close()
         
    
    DATA_Load=sph.ScalarPatternUniform(load_data_reverse[:,0:-1])
    clm=sph.spht(DATA_Load,Degree,Degree)
##LOAD COMPONENT SET TO ZERO FOR GRAVITY PERTURBATION
#        clm[0,0]=0
#        clm[1,:]=[0,0,0]
    
#==============================================================================         
                
    for j in range(N_layer):
        if N_layer==1:
            flag_cis=2
        elif j==0:
            flag_cis=0
        elif j<N_layer-1:
            flag_cis=1
        else:
            flag_cis=2
        
        #Clm  = clm.copy()*0
        Zero=sph.ScalarPatternUniform(np.zeros([grid1+1,2*grid1]))
        Clm  = sph.spht(Zero,Degree,Degree)
        
        if flag_cis==0:
                           
            
            for ii in range(Degree+1):

                fac=4*pi*G/(2*ii+1)*Radius[0]
                Clm[ii,:]=clm[ii,:]*fac
                  
                for k in range(N_layer):
                    fac=4*pi*G/(2*ii+1)*Radius[0]*(Density[k]-Density[k+1])*(Radius[0]/Radius[k])**(ii-1)
                    Clm[ii,:]=Clm[ii,:]+fac*clm_stock[k][ii,:]

                         
        elif flag_cis==1:
            for ii in range(Degree+1):
                fac=4*pi*G/(2*ii+1)*Radius[j]
                Clm[ii,:]=clm[ii,:]*fac
                 
                for k in range(j+1):
                    fac=4*pi*G/(2*ii+1)*Radius[k]*(Density[k]-Density[k+1])*(Radius[k]/Radius[j])**(ii+1)
                    Clm[ii,:]=Clm[ii,:]+fac*clm_stock[k][ii,:]
                 
                for k in range(j+1,N_layer):
                    fac=4*pi*G/(2*ii+1)*Radius[j]*(Density[k]-Density[k+1])*(Radius[j]/Radius[k])**(ii-1)
                    Clm[ii,:]=Clm[ii,:]+fac*clm_stock[k][ii,:]
             
##                
        else: ##Surface layer
            for ii in range(Degree+1):
                fac=4*pi*G/(2*ii+1)*Radius[-1]
                Clm[ii,:]= clm[ii,:]*fac
                                          
                for k in range(N_layer):
                    fac=4*pi*G/(2*ii+1)*Radius[k]*(Density[k]-Density[k+1])*(Radius[k]/Radius[-1])**(ii+1)
                    Clm[ii,:]=Clm[ii,:]+fac*clm_stock[k][ii,:]
    

#==============================================================================!!!MAYBE REVISE FORCES!!
        if DisableSH_0_1==1:
            Clm[0,0]=0
            Clm[1,:]=[0,0,0]

        Phi=np.real(sph.ispht(Clm,grid1+1,grid1*2).cdata)
        Phi=np.hstack([Phi, Phi[:,0:1]])
        
        exec("Data_in"+str(j)+"[(grid1+1)*i:(grid1+1)*(i+1),:]=Phi")
##            if j == (N_layer-1):
##                geoid_normal                             = Phi/Gacc[j]
##                geoid[(grid1+1)*i:(grid1+1)*(i+1),:]     = geoid_normal
    
    CoG_movement_lm_vector = np.zeros([3])

    CoG_movement_lm = np.sqrt(3.0/np.pi)*Clm[1,:]/Gacc[N_layer-1]
    CoG_movement_lm_vector[0] = 1.0/np.sqrt(2.0)*np.imag(CoG_movement_lm[2])
    CoG_movement_lm_vector[1] = 1.0/2.0*np.real(CoG_movement_lm[1])
    CoG_movement_lm_vector[2] = -1.0/np.sqrt(2.0)*np.real(CoG_movement_lm[2])
    
    r_cog_vec_Steps[i,:] = [CoG_movement_lm_vector[0], CoG_movement_lm_vector[1], CoG_movement_lm_vector[2]]
            
############################################################################################################
############################This is outside the iteration procedure and so it comes after it ###############   
#==============================================================================
#         MoI[i,:]=sum(temp_moi)
#         K2[0,i]=np.real(Clm[2,0])
        
#==============================================================================
    eval("np.savetxt(os.path.join(complete_dir, 'CoG_vector.dat'),r_cog_vec_Steps)")
    eval("np.savetxt(os.path.join(complete_dir, 'Check_geoid.dat'),geoid)")


    # Here we save the deflections, so this is really important as saving the stresses to be added to the GIA ones has
    # to happen in a similar fashion.
    for i in range(N_step):
        U = Deflection_Surface[:,:,i]
        Deflection_surface[(i+1)*(grid1+1):(i+2)*(grid1+1),:] = Deflection_Surface[:,:,i]
        eval("np.savetxt(os.path.join(complete_dir, 'Deflection"+str(i)+".dat'),U)")

        dx = dx_surface[:,:,i]
        eval("np.savetxt(os.path.join(complete_dir, 'dx_surface"+str(i)+".dat'),dx)")
                
        dy = dy_surface[:,:,i]
        eval("np.savetxt(os.path.join(complete_dir, 'dy_surface"+str(i)+".dat'),dy)")
                
        dz = dz_surface[:,:,i]
        eval("np.savetxt(os.path.join(complete_dir, 'dz_surface"+str(i)+".dat'),dz)")
        
#         for j in range (N_layer):
#             Deflection_layers[(j)*(grid1+1):(j+1)*(grid1+1),:] = Deflection[:,:,j]
#             if j == (N_layer-1):       
#                 eval("np.savetxt('Deflection_layers_step"+str(i)+".dat',Deflection_layers)")
    eval("np.savetxt(os.path.join(complete_dir, 'Check_deflection.dat'),Deflection_surface)")
   
    for j in range(N_layer):
        eval("np.savetxt(os.path.join(complete_dir, r'Data_in"+str(j)+".dat'),Data_in"+str(j)+")")
        
## delete printing statement
##        U_l = Deflection[:,:,j] 
##        eval("np.savetxt('Deflection_layer"+str(j)+".dat',U_l)")

    if SLE_true == 1 and Only_Eustatic == 0:

        eta_2 = 1E-4
        Convergence_Zeta0 = 1
        
#         Deflection_surface_normal = DataIce*0.0
#         geoid_normal = DataIce*0.0
        
        Delta_ice           = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
        #Deflection_surface  = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
        #geoid               = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
        
        ####################
        Delta_ice[(grid1+1):(N_step+1)*(grid1+1),:]             = DataIce
#         Deflection_surface[(grid1+1):(N_step+1)*(grid1+1),:]    = Deflection_surface_normal
        exec("phi1_matrix = np.loadtxt('Data_in"+str(N_layer-1)+".dat')")
        geoid[(grid1+1):(N_step+1)*(grid1+1),:]                 = phi1_matrix/Gacc[N_layer-1]
        ####################

        print ('Starting SLE')
        DATA=sph.ScalarPatternUniform(Topo_present[:,0:-1])
        Topo_present_lm=sph.spht(DATA,Degree,Degree)
        
        Ice_present = DataIce[N_step*grid1:(N_step+1)*grid1,:]
        C_ocean_present = np.tile(Topo_present, (N_step+1,1)) < 0 
        ## Create initial topo and ocean surface matrices
        Topo       = np.tile(Topo_present, (N_step+1,1))
        
        t_convergence = -1
        eta1=0
        while eta1 < 0.0005:
            t_convergence = t_convergence+1
            DATA_Delta_ice_initial=sph.ScalarPatternUniform(Delta_ice[t_convergence*(grid1+1):(t_convergence+1)*(grid1+1),0:-1])
            Delta_ice_initial_lm = sph.spht(DATA_Delta_ice_initial,Degree,Degree)
            
            DATA_Delta_ice_present=sph.ScalarPatternUniform(Delta_ice[N_step*(grid1+1):(N_step+1)*(grid1+1),0:-1])
            Delta_ice_present_lm = sph.spht(DATA_Delta_ice_present,Degree,Degree)
            
            Delta_ice_initial_sum =1
            Delta_ice_present_sum =1
            for ii in range(Degree+1):
                Delta_ice_initial_sum = Delta_ice_initial_sum + np.sum(np.absolute(Delta_ice_initial_lm[ii,:]))
                Delta_ice_present_sum = Delta_ice_present_sum + np.sum(np.absolute(Delta_ice_present_lm[ii,:]))
            
            eta1 = np.absolute((Delta_ice_initial_sum - Delta_ice_present_sum)/Delta_ice_present_sum)
            print ('t_convergence')
            print (t_convergence)
            print (eta1)
        k=1
        Convergence_Zeta0_preprevious = 100
        Convergence_Zeta0_previous = 10
        Delta_background_geoid_lastK = np.zeros([N_step+1],dtype=complex64)
        Delta_background_geoid_lastK2 = np.zeros([N_step+1],dtype=complex64)

        while Convergence_Zeta0>eta_2: 
        
            j=0
            Delta_SL = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            IceGround_1 = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            IceGround_2 = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            IceGround_3 = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            IceGround_area = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            Ice_Grounded = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            Beta_NoIce = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            C_IceFreeOcean = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            Eustatic_map = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            d_S = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
            
            if (np.absolute(Convergence_Zeta0/Convergence_Zeta0_preprevious-1)<eta_2*3333):
                added_threshold = (Delta_background_geoid_lastK2-Delta_background_geoid_lastK)/2
                for j in range(N_step):
                    Topo[(grid1+1)*j:(grid1+1)*(j+1),:] = Topo[(grid1+1)*j:(grid1+1)*(j+1),:] - np.real(added_threshold[j] - added_threshold[N_step])
                j=0
            C_ocean = Topo < 0
            IceGround_1 = C_ocean == False
            if initial_ice == 1: 
                IceGround_2 = (Delta_ice + np.tile(Initial_Ice, (N_step+1,1))) > (abs(Topo)*Density_water/Density_Ice)
            else:
                IceGround_2 = Delta_ice > 0
            IceGround_3 = IceGround_2*C_ocean
            IceGround_area = (IceGround_1*1+IceGround_3*1)>0
            Ice_Grounded = Delta_ice* IceGround_area
            Beta_NoIce = Ice_Grounded < 0.00001
            C_IceFreeOcean = Beta_NoIce*C_ocean
            Bathemtry_Topo = np.tile(Topo_present, (N_step+1,1))*C_ocean_present*Beta_NoIce
            
            ## Just here to comply with Kendall et al. convention of glacial cycle counts
            C_IceFreeOcean_lastK = C_IceFreeOcean
            Ice_Grounded_lastK = Ice_Grounded
            Delta_SL_spatial = geoid - Deflection_surface## Must be loaded for each point in time
        
            RO = Delta_SL_spatial*C_IceFreeOcean_lastK
            TO = np.tile(Topo[0:(grid1+1),:], (N_step+1,1))*(C_IceFreeOcean*1-np.tile(C_IceFreeOcean[0:(grid1+1),:], (N_step+1,1))*1)
        
            for j in range(N_step+1):      
        #        i=0
                if j==0:
                    DATA_0=sph.ScalarPatternUniform(np.zeros([grid1+1,grid1*2]))
                    Delta_S_lm = sph.spht(DATA_0,Degree,Degree)
                    Delta_S = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
                    Eustatic = np.zeros([(N_step+1)*(grid1+1),grid1*2+1])
                    d_S_lm = list()
                else:
                    
                    Delta_ice_grounded_lastK = Ice_Grounded[(grid1+1)*j:(grid1+1)*(j+1),:]
                    DATA_Delta_ice_grounded=sph.ScalarPatternUniform(Delta_ice_grounded_lastK[:,0:-1])
                    Delta_ice_grounded_lm_lastK = sph.spht(DATA_Delta_ice_grounded,Degree,Degree)
                    
                    DATA_C_iceFree=sph.ScalarPatternUniform(C_IceFreeOcean_lastK[(grid1+1)*j:(grid1+1)*(j+1),0:-1])
                    C_IceFreeOcean_lastK_J_lm = sph.spht(DATA_C_iceFree,Degree,Degree)
                    
                    DATA_TO=sph.ScalarPatternUniform(TO[(grid1+1)*j:(grid1+1)*(j+1),0:-1])
                    TO_lm = sph.spht(DATA_TO,Degree,Degree)
                 
                        
                    DATA_RO_J=sph.ScalarPatternUniform(RO[(grid1+1)*j:(grid1+1)*(j+1),0:-1])
                    RO_J_lm = sph.spht(DATA_RO_J,Degree,Degree)
        
                    Delta_background_geoid_J_lm = 1/C_IceFreeOcean_lastK_J_lm[0,0]*(-Density_Ice/Density_water*Delta_ice_grounded_lm_lastK[0,0] - RO_J_lm[0,0] + TO_lm[0,0] )
                    Delta_background_geoid_lastK2[j] = Delta_background_geoid_lastK[j]
                    Delta_background_geoid_lastK[j] = Delta_background_geoid_J_lm 
                    
                    Eustatic_map[(grid1+1)*j:(grid1+1)*(j+1),:] = np.real(Delta_background_geoid_J_lm)*C_IceFreeOcean[(grid1+1)*j:(grid1+1)*(j+1),:]
            
                    Delta_SL[(grid1+1)*j:(grid1+1)*(j+1),:] = Delta_SL_spatial[(grid1+1)*j:(grid1+1)*(j+1),:] + np.real(Delta_background_geoid_J_lm)#Eustatic_val
            
            if k==1:
                DATA_Topo_initial=sph.ScalarPatternUniform(Topo[t_convergence*(grid1+1):(t_convergence+1)*(grid1+1),0:-1])
                Topo_initial_lm_LastK = sph.spht(DATA_Topo_initial,Degree,Degree)
                
            Topo = np.tile(Topo_present + Delta_SL[(N_step)*(grid1+1):(N_step+1)*(grid1+1),:], (N_step+1,1))-Delta_SL #check line
            
            if initial_ice == 1:
                Delta_S =  RO + Eustatic_map - TO
            else:
                Delta_S =  (RO + Eustatic_map - TO)*(np.ones([(N_step+1)*(grid1+1),grid1*2+1])-IceGround_3)
            
            DATA_Topo_initial=sph.ScalarPatternUniform(Topo[t_convergence*(grid1+1):(t_convergence+1)*(grid1+1),0:-1])
            Topo_initial_lm = sph.spht(DATA_Topo_initial,Degree,Degree)
            
            Topo_initial_lm_sum =0
            Topo_initial_lm_LastK_sum =0
            for ii in range(Degree+1):
                Topo_initial_lm_sum = Topo_initial_lm_sum + np.sum(np.absolute(Topo_initial_lm[ii,:]))
                Topo_initial_lm_LastK_sum =Topo_initial_lm_LastK_sum + np.sum(np.absolute(Topo_initial_lm_LastK[ii,:]))
            
            Convergence_Zeta0_preprevious = Convergence_Zeta0_previous
            Convergence_Zeta0_previous = Convergence_Zeta0
            Convergence_Zeta0 = np.absolute((Topo_initial_lm_sum - Topo_initial_lm_LastK_sum)/Topo_initial_lm_LastK_sum)
            print ('Convergence criterion value:')
            print (Convergence_Zeta0)
            
            Topo_initial_lm_LastK = Topo_initial_lm
            Topo_initial_lm = Topo_initial_lm*0
            
            d_S_lastK = d_S
            k=k+1
           
        np.savetxt(os.path.join(complete_dir, 'Data_in_Slvl.dat'),Delta_S[(grid1+1):(N_step+1)*(grid1+1),:])
        np.savetxt(os.path.join(complete_dir, 'RSL.dat'),Delta_SL)
        np.savetxt(os.path.join(complete_dir, 'IceFreeOcean.dat'),C_IceFreeOcean)
        np.savetxt(os.path.join(complete_dir, 'RSL_spatial.dat'),Delta_SL_spatial)
        
        print ('done with SLE')
    
