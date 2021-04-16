This README will provide a short overview of the scripts used in ABAQUS to perform model creation and stress combination. The two main scripts to run, which call a series of other scripts, are model_gen and iter_ult. The first one is what actually creates the .cae file, the base mode on which to operate. The second one applies forces due to the ice load and other options such as centrifugal ones deriving from the liquid core if the option is selected and calculates stresses and deformations. Its is also here that the stress coupling procedure takes place.

All these scripts have only been modified by me and have been originally written by Bas Blank and Haiyang Hu of the Delft University of Technology.

## Mmodel_data2_top
This first script defines a series of parameters used by all the others, ranging from the density of water and ice to Boolean constants related to loading options.

## Model_data2
This script defines the material constants for each of the solid Earth's layers and assigns them their density. After this, either the ice load is interpolated in a grid using scatterLoad2blockLoad or loads the initial world topography plus the ice sheet and modifies the ice load at the last 2 time steps if the elastic filter option is on.

### ScatterLoad2BlockLoad
This is a subscript of model_data2 and is called from there. Its function is to mesh and then interpolate the ice load data over a grid from a series of initial scattered points. This is done for every ice history time step.

## Model_gen_cone_rampV6_3

## Initial_CoG_correction

## Iter_ult_v7_3SLE

## sph_tools_TPW7_4SLE

## force_definition_routine_v2
This script is not part of the standard scripts given at the beginning of the thesis project but it has been developed as a base for the use of more realistic stresses in the iter_ult script. The stresses are scaled based on the plate velocity maps that can be found in literature and this script is supposed to define an average stress for the AOI from velocity maps to then scale this value with depth for every element.

This script has already been implemented in iter_ult and commented out, but it has also been presented as a separate script in this folder for better modification. 

