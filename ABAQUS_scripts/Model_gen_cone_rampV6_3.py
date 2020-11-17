  ## A multi-layer model generator for Abaqus-CAE
## All parameters are set in Model_data.py
## Author1: Bas Blank
## Author2: Haiyang Hu
## Last Build : 10-04-2018

from abaqus import *
from abaqusConstants import *
from caeModules import *
from odbAccess import*
import numpy as np
import os


Dir=r"/home/fabri/Earth_model_abaqus_SLE0/"
Dir2 = r"/home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2"
sys.path.insert(0,Dir2)
sys.path.insert(1,Dir)
import Model_data2_top
from Model_data2_top import *
import Model_data2
from Model_data2 import *
os.chdir(complete_dir)
print "The folder for Model_gen is: "
print os.getcwd()
#############################################################################################################
#### Model Building 

Mdb()
mdb.models.changeKey('Model-1',Model_name)

##  Generate Core
s = mdb.models[Model_name].ConstrainedSketch(name='__profile__',sheetSize=10000.0)
g, v, d= s.geometry, s.vertices, s.dimensions
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -5000.0), point2=(0.0, 5000.0))
s.FixedConstraint(entity=g[2])

s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, Radius[0]), point2=(0.0, -Radius[0]), direction=CLOCKWISE)

s.Line(point1=(0.0, -Radius[0]), point2=(0.0, Radius[0]))
p = mdb.models[Model_name].Part(name='L0', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models[Model_name].parts['L0']
p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)

p.PartitionCellByPlaneThreePoints(point1=(0,0,0),point2=(0,1e3,0), point3=(0,0,1e3),cells=p.cells)                                                     
p.PartitionCellByPlaneThreePoints(point1=(0,0,0),point2=(0,1e3,0), point3=(1e3,0,0),cells=p.cells)
p.PartitionCellByPlaneThreePoints(point1=(0,0,0),point2=(0,0,1e3), point3=(1e3,0,0),cells=p.cells)   

v=p.vertices

if N_layer==1:
        p.Set(faces=p.faces.findAt((0,Radius[0],0),).getFacesByFaceAngle(10), name='IF0')    #node set of inner layer which is outside of core (Interface0)                                             
p.Set(vertices=v.findAt(((0,0,0),)),name='Center')   #Center of core                               
p.Surface(side1Faces=p.faces.findAt((0,Radius[0],0),).getFacesByFaceAngle(10), name='Sout0') #Core Surface                          
                              
                              
##  Generate Outer Layers
for i in range(1,N_layer):                                
    s = mdb.models[Model_name].ConstrainedSketch(name='__profile__',    sheetSize=10000.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.ConstructionLine(point1=(0.0, -5000.0), point2=(0.0, 5000.0))
 
    s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, Radius[i-1]), point2=(0.0, -Radius[i-1]), direction=CLOCKWISE)
    s.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, Radius[i]), point2=(0.0, -Radius[i]), direction=CLOCKWISE)
    s.Line(point1=(0.0, -Radius[i-1]), point2=(0.0, -Radius[i]))                                                    
    s.Line(point1=(0.0, Radius[i-1]), point2=(0.0, Radius[i]))
    
    p = mdb.models[Model_name].Part(name='L'+str(i), dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
    #s.Line(point1=(Radius[i-1]*sin(alpha), -Radius[i-1]*cos(alpha)), point2=(Radius[i]*sin(alpha), -Radius[i]*cos(alpha)))
    s.unsetPrimaryObject()
    del mdb.models[Model_name].sketches['__profile__']
                                                    
    p.PartitionCellByPlaneThreePoints(point1=(0,0,0),point2=(0,1e3,0), point3=(0,0,1e3),cells=p.cells)                                                     
    p.PartitionCellByPlaneThreePoints(point1=(0,0,0),point2=(0,1e3,0), point3=(1e3,0,0),cells=p.cells)
    p.PartitionCellByPlaneThreePoints(point1=(0,0,0),point2=(0,0,1e3), point3=(1e3,0,0),cells=p.cells)   
    
    if i==1:
        p.Set(faces=p.faces.findAt((0,Radius[0],0),).getFacesByFaceAngle(10), name='IF0')    #node set of inner layer which is outside of core (Interface0)   
    p.Set(faces=p.faces.findAt((0,Radius[i],0),).getFacesByFaceAngle(10), name='IF'+str(i))    #node set of outer layer                                                   
    p.Surface(side1Faces=p.faces.findAt((0,Radius[i-1],0),).getFacesByFaceAngle(10), name='Sin'+str(i)) #Inner L1 Surface                                                    
    p.Surface(side1Faces=p.faces.findAt((0,Radius[i],0),).getFacesByFaceAngle(10), name='Sout'+str(i)) #Outter L1 Surface          
   
##Material, section

for i in range(N_layer):
    mdb.models[Model_name].Material(name='M'+str(i))
    mdb.models[Model_name].materials['M'+str(i)].Elastic(moduli=INSTANTANEOUS,table=((Y_mod[i], Poi), ))
    mdb.models[Model_name].materials['M'+str(i)].Density(table=((Density[i], ), ))
    if Vis[i]!=0 and Vis[i]<1e40:
##        mdb.models[Model_name].materials['M'+str(i)].Creep(law=STRAIN,table=((1.0/3/Vis[i], 1, 0.0),))
        mdb.models[Model_name].materials['M'+str(i)].Viscoelastic(domain=TIME, time=PRONY,     table=((1-1e-10, 0.0, Vis[i]/Y_mod[i]*3), ))


    mdb.models[Model_name].HomogeneousSolidSection(name='S'+str(i),  material='M'+str(i), thickness=None)
    p = mdb.models[Model_name].parts['L'+str(i)]   
    p.SectionAssignment(region=regionToolset.Region(cells=p.cells), sectionName='S'+str(i), offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)   

    
## Assemblely and mesh
a = mdb.models[Model_name].rootAssembly
if flag_lc==0:
    start=0
else:
    start=1
        
    a.DatumCsysByThreePoints(name='Coor', coordSysType=SPHERICAL, origin=(0.0, 0.0, 0.0), line1=(0.0, 0.0, 1.0), line2=(1.0, 0.0, 0.0))

if N_layer>2 or (N_layer==2 and flag_lc==0):
    for i in range(0,N_layer):
        p = mdb.models[Model_name].parts['L'+str(i)]
        a.Instance(name='I'+str(i), part=p, dependent=ON)
        if i==0:
            region = a.instances['I0'].sets['Center']
    if Coarse_Layers == 1:
##        p = mdb.models[Model_name].parts['L1']
##        a.Instance(name=Model_name + '_Low', part=p, dependent=ON)
        mdb.models[Model_name].parts.changeKey('L1',Model_name + '_Low')
        p = mdb.models[Model_name].parts[Model_name + '_Low']
        a.Instance(name=Model_name + '_Low', part=p, dependent=ON)
    else:       
        a.InstanceFromBooleanMerge(name=Model_name + '_Low', instances=list(a.instances['I'+str(j)] for j in range(1,Coarse_Layers+1)), keepIntersections=ON, originalInstances=SUPPRESS, mergeNodes=ALL, nodeMergingTolerance=1e-06, domain=BOTH)
        a.instances.changeKey(Model_name+'_Low-1',Model_name + '_Low')
    a.InstanceFromBooleanMerge(name=Model_name, instances=list(a.instances['I'+str(j)] for j in range(Coarse_Layers+1,N_layer)), keepIntersections=ON, originalInstances=SUPPRESS, mergeNodes=ALL, nodeMergingTolerance=1e-06, domain=BOTH)
    a.instances.changeKey(Model_name+'-1',Model_name)
elif N_layer==2:
    mdb.models[Model_name].parts.changeKey('L1',Model_name)
    p = mdb.models[Model_name].parts[Model_name]
    a.Instance(name=Model_name, part=p, dependent=ON)
else:
    mdb.models[Model_name].parts.changeKey('L0',Model_name)
    p = mdb.models[Model_name].parts[Model_name]
    a.Instance(name=Model_name, part=p, dependent=ON)

elemType1 = mesh.ElemType(elemCode=C3D8RH, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D6H, elemLibrary=STANDARD,     secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType3 = mesh.ElemType(elemCode=C3D4H, elemLibrary=STANDARD)


p = mdb.models[Model_name].parts[Model_name]
c = p.cells
p.setMeshControls(regions=c,elemShape=HEX_DOMINATED, technique=SWEEP,  algorithm=ADVANCING_FRONT)


#### edge 1

##Top of surface mesh adjustment
p.seedPart(size=Seeds, deviationFactor=0.05, minSizeFactor=0.1)

Par = PolarCirkelRadiusDeg/90.0
#p = mdb.models[Model_name].parts[Model_name]
e = p.edges

pickedEdges1=e.getClosest(coordinates=((Radius[-1]*np.sin((PolarCirkelRadiusDeg/2)/180*np.pi),-Radius[-1]*np.cos((PolarCirkelRadiusDeg/2)/180*np.pi),0),))[0][0]
p.PartitionEdgeByParam(edges=pickedEdges1, parameter=Par)

for i in range(Coarse_Layers+1,N_layer):
    if i==Coarse_Layers+1:
        pickedEdges=e.getClosest(coordinates=((Radius[1]*np.sin((PolarCirkelRadiusDeg/2)/180*np.pi),-Radius[Coarse_Layers]*np.cos((PolarCirkelRadiusDeg/2)/180*np.pi),0),))[0][0]
        p.PartitionEdgeByParam(edges=pickedEdges, parameter=Par)

    pickedEdges=e.getClosest(coordinates=((Radius[i]*np.sin(alpha),-Radius[i]*np.cos(alpha),0.0),))[0][0]
    
    if i<(N_layer-1):
        p.PartitionEdgeByParam(edges=pickedEdges, parameter=Par)
    
    R_av = (Radius[i-1]+Radius[i])/2
    v1= p.vertices.getClosest(coordinates=((Radius[i-1]*sin(alpha), -Radius[i-1]*cos(alpha),0,),))[0][0]
    v2= p.vertices.getClosest(coordinates=((Radius[i]*sin(alpha), -Radius[i]*cos(alpha),0,),))[0][0]
    F1 = p.faces.getClosest(coordinates =((R_av*np.sin(alpha),-R_av*np.cos(alpha),0.0,),))[0][0]
    p.PartitionFaceByShortestPath(point1=v1,point2=v2, faces=F1)
    

    pickedCells = p.cells.getByBoundingSphere((0.0,0.0,0.0), 1.1*Radius[-1])
    pickedEdges = e.getClosest(coordinates=((R_av*np.sin(alpha),-R_av*np.cos(alpha),0.0,),))[0][0]
    sweepEdge = e.findAt((Radius[i-1]*np.cos(1.0/4.0*np.pi),0.0,-Radius[i-1]*np.sin(1.0/4.0*np.pi)))
    p.PartitionCellBySweepEdge(sweepPath=sweepEdge, cells=pickedCells, edges=(pickedEdges, ))
    
    pickedCells = p.cells.getByBoundingSphere((0.0,0.0,0.0), 1.1*Radius[-1])
    pickedEdges = e.getClosest(coordinates=((0.0,-R_av*np.cos(alpha),-R_av*np.sin(alpha),),))[0][0]
    sweepEdge = e.findAt((-Radius[i-1]*np.cos(1.0/4.0*np.pi),0.0,-Radius[i-1]*np.sin(1.0/4.0*np.pi)))
    p.PartitionCellBySweepEdge(sweepPath=sweepEdge, cells=pickedCells, edges=(pickedEdges, ))
    
    pickedCells = p.cells.getByBoundingSphere((0.0,0.0,0.0), 1.1*Radius[-1])
    pickedEdges = e.getClosest(coordinates=((-R_av*np.sin(alpha),-R_av*np.cos(alpha),0.0,),))[0][0]
    sweepEdge = e.findAt((-Radius[i-1]*np.cos(1.0/4.0*np.pi),0.0,Radius[i-1]*np.sin(1.0/4.0*np.pi)))
    p.PartitionCellBySweepEdge(sweepPath=sweepEdge, cells=pickedCells, edges=(pickedEdges, ))
    
    pickedCells = p.cells.getByBoundingSphere((0.0,0.0,0.0), 1.1*Radius[-1])
    pickedEdges = e.getClosest(coordinates=((0.0,-R_av*np.cos(alpha),R_av*np.sin(alpha),),))[0][0]
    sweepEdge = e.findAt((Radius[i-1]*np.cos(1.0/4.0*np.pi),0.0,Radius[i-1]*np.sin(1.0/4.0*np.pi)))
    p.PartitionCellBySweepEdge(sweepPath=sweepEdge, cells=pickedCells, edges=(pickedEdges, ))
    
    pickedEdges = e.findAt((0.0,-R_av,0.0))
    p.seedEdgeBySize(edges=(pickedEdges,), size=R_Target_seed, deviationFactor=0.05,minSizeFactor=0.1, constraint=FINER)
    pickedEdges1 = e.findAt((-Radius[i-1]*np.sin(alpha/2),-Radius[i-1]*np.cos(alpha/2),0))
    pickedEdges2 = e.findAt((Radius[i-1]*np.sin(alpha/2),-Radius[i-1]*np.cos(alpha/2),0))
    pickedEdges3 = e.findAt((0,-Radius[i-1]*np.cos(alpha/2),-Radius[i-1]*np.sin(alpha/2)))
    pickedEdges4 = e.findAt((0,-Radius[i-1]*np.cos(alpha/2),Radius[i-1]*np.sin(alpha/2)))
    p.seedEdgeBySize(edges=(pickedEdges1,pickedEdges2,pickedEdges3,pickedEdges4,), size=Plane_Target_seed, deviationFactor=0.05,minSizeFactor=0.1, constraint=FINER)

pickedEdges1 = e.findAt((-Radius[-1]*np.sin(alpha/2),-Radius[-1]*np.cos(alpha/2),0))
pickedEdges2 = e.findAt((Radius[-1]*np.sin(alpha/2),-Radius[-1]*np.cos(alpha/2),0))
pickedEdges3 = e.findAt((0,-Radius[-1]*np.cos(alpha/2),-Radius[-1]*np.sin(alpha/2)))
pickedEdges4 = e.findAt((0,-Radius[-1]*np.cos(alpha/2),Radius[-1]*np.sin(alpha/2)))
p.seedEdgeBySize(edges=(pickedEdges1,pickedEdges2,pickedEdges3,pickedEdges4,), size=Plane_Target_seed, deviationFactor=0.05,minSizeFactor=0.1, constraint=FINER)

pickedEdges = e.getByBoundingCylinder((0.0,-1.1*Radius[-1],0.0),(0.0,-0.9*Radius[0],0.0), 1000.0)
p.seedEdgeBySize(edges=pickedEdges, size=75000.0, deviationFactor=0.05,minSizeFactor=0.1, constraint=FINER)
   
cells=c.getByBoundingBox(-Radius[-1]-100,-Radius[-1]-100,-Radius[-1]-100,Radius[-1]+100,Radius[-1]+100,Radius[-1]+100)
pickedRegions =(cells, )


# generate mesh
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))  
p.generateMesh()
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))  

## generate core
p = mdb.models[Model_name].parts['L0']
c = p.cells
cells=c.getByBoundingBox(-Radius[0]-100,-Radius[0]-100,-Radius[0]-100,Radius[0]+100,Radius[0]+100,Radius[0]+100)
pickedRegions =(cells, )
p.setMeshControls(regions=c,elemShape=HEX, technique=STRUCTURED,  algorithm=DEFAULT)
p.seedPart(size=Seeds*2, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))

#a.Instance(name='Core', part=p, dependent=ON)

a = mdb.models[Model_name].rootAssembly
region1=a.instances['I0'].surfaces['Sout0']
region2=a.instances[Model_name + '_Low'].surfaces['Sin1']
mdb.models[Model_name].Tie(name='Constraint-1', master=region1, slave=region2, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
##session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
##    predefinedFields=ON, interactions=OFF, constraints=OFF, 
##    engineeringFeatures=OFF)
#a = mdb.models[Model_name].rootAssembly
#region = a.instances['Core'].sets['Center']

# for lower mantle

a = mdb.models[Model_name].rootAssembly
region1=a.instances[Model_name+'_Low'].surfaces['Sout'+ str(Coarse_Layers)] 
a = mdb.models[Model_name].rootAssembly
region2=a.instances[Model_name].surfaces['Sin' + str(Coarse_Layers+1)]
mdb.models[Model_name].Tie(name='Constraint-2', master=region1, slave=region2, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
p = mdb.models[Model_name].parts[Model_name + '_Low']
p.setMeshControls(regions=c,elemShape=HEX, technique=STRUCTURED,  algorithm=DEFAULT)
p.seedPart(size=Seeds*LM_res, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
p = mdb.models[Model_name].parts[Model_name]

if fixed_point ==1:
    ## Boundary condition.
    a = mdb.models[Model_name].rootAssembly
    region_core = a.instances['I0'].sets['Center']
    mdb.models[Model_name].EncastreBC(name='BC-1', createStepName='Initial',     region=region_core, localCsys=None) ##  fixing center of Earth

if hires_rotate_grid == 1:
    a.rotate(instanceList=('Earth_Low', 'Earth'), axisPoint=(0.0, 0.0, 0.0), axisDirection=(-1, 0.0, 0.0), angle=(phi_hr+90))
    a.rotate(instanceList=('Earth_Low', 'Earth'), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1, 0.0), angle=theta_hr) 
                                                          
## Steps                       
for i in range(N_step):                         
    if i==0:
	if Rampload_enabled==1:
	    mdb.models[Model_name].ViscoStep(name='Step'+str(i+1), previous='Initial',    timePeriod=Time[i], maxNumInc=500, initialInc=Time[i]/1e2,     minInc=Time[i]/1e6, maxInc=Time[i]/25, cetol=CETOL, amplitude=RAMP)
	else:
	    mdb.models[Model_name].ViscoStep(name='Step'+str(i+1), previous='Initial',    timePeriod=Time[i], maxNumInc=500, initialInc=Time[i]/1e2,     minInc=Time[i]/1e6, maxInc=Time[i]/25, cetol=CETOL)
		
    else:
        Gap=Time[i]-Time[i-1]
	if Rampload_enabled==1:
	    mdb.models[Model_name].ViscoStep(name='Step'+str(i+1), previous='Step'+str(i),    timePeriod=Gap, maxNumInc=500, initialInc=Gap/1e2,     minInc=Gap/1e6, maxInc=Gap/25, cetol=CETOL, amplitude=RAMP)
	else:
	    mdb.models[Model_name].ViscoStep(name='Step'+str(i+1), previous='Step'+str(i),    timePeriod=Gap, maxNumInc=500, initialInc=Gap/1e2,     minInc=Gap/1e6, maxInc=Gap/25, cetol=CETOL)

## Foundation
for i in range(N_layer-1):
    if Density[i]>Density[i+1]:
        if i < Coarse_Layers:
            region=a.instances[Model_name + '_Low'].surfaces['Sin'+str(i+1)]
        else:
            region=a.instances[Model_name].surfaces['Sin'+str(i+1)]
        mdb.models[Model_name].ElasticFoundation(name='F'+str(i), createStepName='Initial', surface=region, stiffness=(Density[i]-Density[i+1])*Gacc[i])
##        stiff = (Density[i]-Density[i+1])*Gacc[i]
region=a.instances[Model_name].surfaces['Sout'+str(N_layer-1)]
mdb.models[Model_name].ElasticFoundation(name='F'+str(N_layer-1), createStepName='Initial', surface=region, stiffness=Density[N_layer-1]*Gacc[N_layer-1])

## Many of the similarities with Maaike's scripts end here. Her viscosteps are deinfed in the TPW script.

##stiff = Density[N_layer-1]*Gacc[N_layer-1]

## Define Loads
amp=1.0/3*Omega**2*Radius[-1]**2


if Centrifugal_on < 0 or Centrifugal_on > 1:
    sys.exit("Invalid value for variable P20_on")
    
if Centrifugal_on != 0:
    Coor_sph=mdb.models[Model_name].rootAssembly.datums.values()[0]                ##Define load potential P20(cos(T))
    mdb.models[Model_name].ExpressionField(name='P20', localCsys=Coor_sph, description='', expression='0.5*(3*cos(P)**2-1)')                               

    if N_layer>1:
        region = a.instances[Model_name].surfaces['Sin1']
        mdb.models[Model_name].Pressure(name='Centrifugal0', createStepName='Step1',region=region, distributionType=FIELD, field='P20', magnitude=(-1)**flag_lc*(Density[0]-Density[1])*amp*Radius[0]**2/Radius[-1]**2, amplitude=UNSET)  ##Take care of the sign, *Radius[0]**2/Radius[-1]**2
        for i in range(1,N_layer):
            if Density[i]>Density[i+1]:
                if i<N_layer-1:                
                    region = a.instances[Model_name].surfaces['Sout'+str(i)]                                           
                    mdb.models[Model_name].Pressure(name='Centrifugal'+str(i), createStepName='Step1',region=region, distributionType=FIELD, field='P20', magnitude=(Density[i]-Density[i+1])*amp*Radius[i]**2/Radius[-1]**2, amplitude=UNSET)  ##*Radius[i]**2/Radius[-1]**2
                else:
                    region = a.instances[Model_name].surfaces['Sout'+str(i)]                                           
                    mdb.models[Model_name].Pressure(name='Centrifugal'+str(i), createStepName='Step1',region=region, distributionType=FIELD, field='P20', magnitude=(Density[i]-Density[i+1])*amp*Radius[i]**2/Radius[-1]**2, amplitude=UNSET)  ##*Radius[i]**2/Radius[-1]**2
                    
    else:
        region = a.instances[Model_name].surfaces['Sout0']
        mdb.models[Model_name].Pressure(name='Centrifugal0', createStepName='Step1',region=region, distributionType=FIELD, field='P20', magnitude=(Density[0])*amp, amplitude=UNSET)
        

## Ice Load + Sea Load

#DataIce=np.loadtxt('DataIce_in.dat')
#gridIce=DataIce.shape[1]/2
S_0 = np.zeros(((grid1+1)*N_step,grid1*2+1))
for i in range(N_step):
          
    #### Create iteration-step 0 ice field and Mass Loss
    if No_IceLoss == 0 & SLE_true==1:
        
        ## Create iteration-step 0 Sea level
    
        R_E = Radius[N_layer-1]
        Surface_Area = 4*pi*R_E**2
        
        IceChange = DataIce[(grid1+1)*i:(grid1+1)*(i+1),:]
        stepsize = 180.0/(grid1)
        ## re-evalute gridIce-1 here!!!!
        colats = np.linspace(0,180,grid1+1)
        
        Sphere_area_matrix = np.zeros((grid1+1, grid1*2+1))
        for int in range(grid1+1):
            if int == 0:
                phi_upper = 0.0
            else:
                phi_upper = (colats[int-1] + colats[int])/2*(pi/180)#Deg2Rad
            if int == max(range(grid1+1)):
                phi_lower = pi
            else:
                phi_lower = (colats[int]+colats[int+1])/2*(pi/180)#Deg2Rad
                
            Sphere_area_matrix[int,:] = abs(0.5*Surface_Area*(np.cos(phi_upper)-np.cos(phi_lower))/(grid1*2))

        M_IceLoss = np.sum(np.multiply(IceChange[0:-1],Sphere_area_matrix[0:-1])*Density_Ice)

        Ocean_function = Ocean_Function[(grid1+1)*i:(grid1+1)*(i+1),:]
        Ocean_surface_fraction = np.sum(Ocean_function[0:-1]*Sphere_area_matrix[0:-1])/Surface_Area
                  
        S_0_step = -M_IceLoss/(Density_water*Surface_Area*Ocean_surface_fraction)*Ocean_function
        S_0[(grid1+1)*i:(grid1+1)*(i+1),:] = S_0_step      
        
    if i == (N_step-1):
        ## print sea_level files
        fileout=open('S0.dat','w')
        np.savetxt(fileout,S_0)
        fileout.close()
            
        fileout=open('Data_in_Slvl.dat','w')
        np.savetxt(fileout,S_0)
        fileout.close()
    
    if No_IceLoss  == 0:
        Sea_lvl_allTimes_ice_eq = S_0*Density_water/Density_Ice
    else:
        Sea_lvl_allTimes_ice_eq = 0
        
    DataLoad =  DataIce  +  Sea_lvl_allTimes_ice_eq
     
## Create initial forces in ABAQUS model
if DisableSH_0_1 ==1:
    print 'Start initial sub-Python process to remove SHdegree 0 and 1\n'
    CHECK=os.system('python2.7 /home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2/sph_tools_initial_load_v2.py') ##1
    print 'Finish initial sub-Python process\n'
    
    DataLoad=np.loadtxt('DataLoad_in.dat')
    
elif Compensate_cog ==1:   
    
    print 'Start initial correction for cog movement\n'
    CHECK=os.system('python2.7 /home/fabri/Earth_model_abaqus_SLE0/parallel_run_scripts_2/Initial_CoG_correction.py') ##1
    print 'Finish initial sub-Python process\n'
    
    theta_mat           = np.tile(np.linspace(0,2*np.pi,(2*grid1+1)) ,[grid1+1, 1]) 
    phi_mat             = np.tile(np.transpose([np.linspace(np.pi/2.0,-np.pi/2.0,(grid1+1))]) ,[1, 2*grid1+1])
    
    r_cog_vec_steps     = np.loadtxt('CoG_vector.dat')
    
    cos_gamma_field = np.tile(np.zeros([grid1+1,2*grid1+1]),[N_step,1])

for i in range(N_step):
    
    Coor_sph=mdb.models[Model_name].rootAssembly.datums.values()[0]

    theta=np.linspace(0,180,grid1+1)
    theta=np.hstack((0,theta))
    theta=theta.reshape(theta.size,1)
        
    lam=np.linspace(0,360,grid1*2+1)
    lam=np.hstack((lam[lam<=180],np.array([J-360 for J in lam if J>180])))   #in Abaqus, Th start from 0~pi then -pi~0.
    
    if Compensate_cog ==1 and DisableSH_0_1 ==0 :
        #r_cog_vec = [r_cog_lm[0], r_cog_lm[1], r_cog_lm[2]]
        r_cog_vec   = r_cog_vec_steps[i,:]
        r_cog_norm  = np.linalg.norm(r_cog_vec)
    
        theta_cog   = np.arctan2(r_cog_vec[0],r_cog_vec[2])
        Delta_theta = theta_cog - theta_mat    
        
        #cos_gamma_field[(grid1+1)*i:(grid1+1)*(i+1),:] = r_cog_vec[1]/r_cog_norm*np.sin(phi_mat)+np.sqrt(r_cog_vec[0]**2+r_cog_vec[2]**2)/r_cog_norm*np.cos(phi_mat)*np.cos(Delta_theta)
        if r_cog_norm==0:
            cos_gamma_field_singleStep = np.zeros([grid1+1,grid1*2+1])
        else:    
            cos_gamma_field_singleStep = r_cog_vec[1]/r_cog_norm*np.sin(phi_mat)+np.sqrt(r_cog_vec[0]**2+r_cog_vec[2]**2)/r_cog_norm*np.cos(phi_mat)*np.cos(Delta_theta)
            

        CoG_correction = r_cog_norm*cos_gamma_field_singleStep
        
        CoG_correction = np.vstack((lam,CoG_correction))
        CoG_correction = np.hstack((theta,CoG_correction))
        CoG_correction = tuple([tuple(J) for J in CoG_correction])
        
        mdb.models[Model_name].MappedField(name='CoG_correction_surface' +str(i)+'_field', description='',
           regionType=POINT, partLevelData=False, localCsys=Coor_sph,
           pointDataFormat=GRID, fieldDataType=SCALAR, gridPointPlane=PLANE32,
           gridPointData={str(Radius[N_layer-1]):CoG_correction})
                                           
        mdb.models[Model_name].MappedField(name='CoG_correction_CMB' +str(i)+'_field', description='',
           regionType=POINT, partLevelData=False, localCsys=Coor_sph,
           pointDataFormat=GRID, fieldDataType=SCALAR, gridPointPlane=PLANE32,
           gridPointData={str(Radius[0]):CoG_correction})
                                           
        a = mdb.models[Model_name].rootAssembly
        region = a.instances[Model_name].surfaces['Sout'+str(N_layer-1)]
        mdb.models[Model_name].Pressure(name='CoG_correction_surface' + str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='CoG_correction_surface' +str(i)+'_field', magnitude=Density[N_layer-1]*Gacc[N_layer-1],amplitude=UNSET)    
        
        region = a.instances[Model_name + '_Low'].surfaces['Sin'+str(1)]
        mdb.models[Model_name].Pressure(name='CoG_correction_CMB' + str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='CoG_correction_CMB' +str(i)+'_field', magnitude=(Density[0]-Density[1])*Gacc[0],amplitude=UNSET)          
        
        if i==0:
            region = a.instances[Model_name].surfaces['Sout'+str(N_layer-1)]
            mdb.models[Model_name].Pressure(name='CoG_correction_surface' + str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='CoG_correction_surface' +str(i)+'_field', magnitude=Density[N_layer-1]*Gacc[N_layer-1],amplitude=UNSET)
            region = a.instances[Model_name + '_Low'].surfaces['Sin'+str(1)]
            mdb.models[Model_name].Pressure(name='CoG_correction_CMB' + str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='CoG_correction_CMB' +str(i)+'_field', magnitude=(Density[0]-Density[1])*Gacc[0],amplitude=UNSET)
        else:
            mdb.models[Model_name].TabularAmplitude(name='Amp-'+str(i), timeSpan=TOTAL, smooth=SOLVER_DEFAULT, data=((
               Time[i-1], 1.0), (Time[i], 0.0)))
                
            region = a.instances[Model_name].surfaces['Sout'+str(N_layer-1)]
            mdb.models[Model_name].Pressure(name='CoG_correction_surface' + str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='CoG_correction_surface' +str(i)+'_field', magnitude=Density[N_layer-1]*Gacc[N_layer-1],amplitude=UNSET)    
            
            region = a.instances[Model_name + '_Low'].surfaces['Sin'+str(1)]
            mdb.models[Model_name].Pressure(name='CoG_correction_CMB' + str(i), createStepName='Step'+str(i+1),region=region, distributionType=FIELD, field='CoG_correction_CMB' +str(i)+'_field', magnitude=(Density[0]-Density[1])*Gacc[0],amplitude=UNSET)          
            
        if i>0:
            mdb.models[Model_name].loads['CoG_correction_surface' + str(i-1)].deactivate('Step'+str(i+1))                                    
            mdb.models[Model_name].loads['CoG_correction_CMB' + str(i-1)].deactivate('Step'+str(i+1))
            
                
##    if DisableSH_0_1 ==1:

    if i>0:
        if filter_elastic==1 and i==(N_step-1):
            New_Force_field=np.vstack((lam,DataLoad[(grid1+1)*(i-1):(grid1+1)*i,:]))
        else:
            New_Force_field=np.vstack((lam,DataLoad[(grid1+1)*i:(grid1+1)*(i+1),:]))
        New_Force_field=np.hstack((theta,New_Force_field)) 
        
        Old_Force_field = New_Force_field
        Force_field= tuple([tuple(J) for J in New_Force_field])
        

    else:
        Force_field=np.vstack((lam,DataLoad[(grid1+1)*i:(grid1+1)*(i+1),:]))
        Force_field=np.hstack((theta,Force_field))
        Old_Force_field = Force_field
        Force_field=tuple([tuple(J) for J in Force_field])
        
    ## Create iteration-step 0 ice load
    a = mdb.models[Model_name].rootAssembly
    region = a.instances[Model_name].surfaces['Sout'+str(N_layer-1)]

    if scatterGrid_load==0: 
        mdb.models[Model_name].MappedField(name='TotalLoad' +str(i)+'_field', description='',
           regionType=POINT, partLevelData=False, localCsys=Coor_sph,
           pointDataFormat=GRID, fieldDataType=SCALAR, gridPointPlane=PLANE32,
           gridPointData={str(Radius[N_layer-1]):Force_field})
        mdb.models[Model_name].Pressure(name='TotalLoad_'+str(i), createStepName='Step'+str(i+1),region=region,
           distributionType=FIELD, field='TotalLoad'+str(i)+'_field', magnitude=Density_Ice*Gacc[N_layer-1],amplitude=UNSET)
    else:
        #mdb.models[Model_name].Pressure(name='TotalLoad_'+str(i), createStepName='Step'+str(i+1),region=region, distributionType=USER_DEFINED, magnitude=1.0, amplitude=UNSET)
        
        iceLoadName = 'IceLoad_' + Model_name + str(i+1) +'_field'
        fileIn = IceDir+'ice_load_xyzf_'+str(i+1)+'.txt'
        DataIce_scat=np.loadtxt(fileIn, delimiter=',')
        mdb.models[Model_name].MappedField(name=iceLoadName, description='', regionType=POINT, partLevelData=False,
           localCsys=None, pointDataFormat=XYZ, fieldDataType=SCALAR, xyzPointData=DataIce_scat)
        mdb.models[Model_name].Pressure(name='IceLoad_' + Model_name + str(i+1), createStepName='Step'+str(i+1),
           region=region, distributionType=FIELD, field=iceLoadName, magnitude=Gacc[N_layer-1], amplitude=UNSET)

    if i>0:
        if scatterGrid_load==0:
            mdb.models[Model_name].loads['TotalLoad_'+str(i-1)].deactivate('Step'+str(i+1))
        else:
            mdb.models[Model_name].loads['IceLoad_'+ Model_name +str(i)].deactivate('Step'+str(i+1))
        
## Define Job

##mdb.Job(name=Model_name, model=Model_name, description='', 
##type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
##memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
##explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
##modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine=Dir + 'user_2.f', 
##scratch='', multiprocessingMode=DEFAULT, numCpus=CPUs, numDomains=CPUs, numGPUs=0)

## Define Output
mdb.models[Model_name].fieldOutputRequests['F-Output-1'].setValues(variables=('U','S','MISES','COORD'))
## Restart request similar to what maaike did:
##for i in range(N_step):
##    mdb.models[Model_name].steps['Step'+str(i)].Restart(frequency=300, numberIntervals=0, overlay=ON, timeMarks=OFF) # restart request
##    if i>0:
##        mdb.models[Model_name].setValues(restartJob=Model_name+str(i-1), restartStep='Step'+str(i-1))    
##    mdb.models[Model_name].fieldOutputRequests['F-Output-1'].setValues(variables=('U','S','MISES','COORD')) # displacement output
##    if len(mdb.models[Model_name].historyOutputRequests)>0:
##        del mdb.models[Model_name].historyOutputRequests['H-Output-1']
# The jobs must then be defined in Iter_ult.
## Save model database
mdb.saveAs(pathName=Model_name+'.cae')
