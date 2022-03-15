import os
import sys
import time
import numpy as np
import random
import matplotlib.pyplot as plt

import irs
import abm

from constants import *
import input
import plots
import cosmology as cosmo
import helper_functions as hf


#### Setting up the random seed
random.seed(input.user_seed)  	
np.random.seed(input.user_seed)  	


#### Observer, Lens and Source redshifts
zO = 0.0
zL = input.zol
zS = input.zos


#### Relevant Angular Diameter Distances (ADD)
Dol = cosmo.angDist(zO, zL)
Dls = cosmo.angDist(zL, zS)
Dos = cosmo.angDist(zO, zS)


#### Critical density (in MASS_SUN/arcsec^2)
critDen = hf.critical_density(Dol,Dls,Dos)


#### Stellar density corresponding to kappa_micro (in MASS_SUN/pc^2)
stelDen = hf.stellar_density(Dol, critDen)
if stelDen < 1e-3:
   sys.exit('Stellar density is very low. \n Please check the parameters.')


#### Getting Einstein angle for one Solar mass lens (in arcseconds)
tEin_scale = hf.Einstein_Angle(Dol,Dls,Dos,input.M_scale)


#### Calculating the relevent source plane region (in arcseconds)
source_plane_FOV = hf.source_plane_FOV(Dos)


#### Calculating the relevant image plane region (in arcseconds)
image_plane_FOV = hf.image_plane_FOV(source_plane_FOV)


#### Generating the microlens population (position (in arcsec) and masses (in solar mass))
ml_position, ml_mass = hf.generate_microlens_population(image_plane_FOV, critDen)


#### Scaling microlens positions (in units of tEin_scale)
ml_position = ml_position/tEin_scale


#### Starting the simulation
#### If IRS is chosen as preferred method in input.py file
if input.method == 'IRS':
    irs.irs_main(Dol, Dls, Dos, tEin_scale, source_plane_FOV, ml_position, ml_mass)
    
#### If ABM is chosen as preferred method in input.py file
if input.method == 'ABM':
    abm.abm_main(Dol, Dls, Dos, tEin_scale, source_plane_FOV, ml_position, ml_mass)


#### Writing the relevant input parametrs
plots.parameters(ml_position, ml_mass, tEin_scale, stelDen)


#### End of the program
