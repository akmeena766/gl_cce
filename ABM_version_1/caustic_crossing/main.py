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
import macrolens as sl
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


#### Stellar density (in MASS_SUN/pc^2) and corresponding kappa_micro
if input.stellar_density < 1e-3:
   sys.exit('Stellar density is very low. \n Please check the parameters.')

kappa_micro = hf.stelDen_to_kappaMicro(input.stellar_density, Dol, critDen)


#### Getting Einstein angle for one Solar mass lens (in arcseconds)
tEin_scale = hf.Einstein_Angle(Dol,Dls,Dos,input.M_scale)


#### Calculating the relevent source plane region (in arcseconds)
source_plane_FOV = hf.source_plane_FOV(Dos)


#### Calculating the relevant image plane region (in arcseconds)
image_plane_FOV = [input.xbox, input.ybox]  


#### Rotation angle of image plane to align arc along y-axis
#kv, g1v, g2v = sl.macrolens_kgamma_NSIE(Dol, Dls, Dos, 0, kappa_micro, 0, 0)
arc_angle = 0.


#### Generating the microlens population (position (in arcsec) and masses (in solar mass))
ml_position, ml_mass = hf.generate_microlens_population(kappa_micro, image_plane_FOV, critDen, arc_angle)


#### Scaling microlens positions (in units of tEin_scale)
ml_position = ml_position/tEin_scale


#### Starting the simulation
#### If IRS is chosen as preferred method in input.py file
if input.method == 'IRS':
    irs.irs_main(Dol, Dls, Dos, tEin_scale, source_plane_FOV, arc_angle, kappa_micro, ml_position, ml_mass)
    
    
#### If ABM is chosen as preferred method in input.py file
if input.method == 'ABM':
    abm.abm_main(Dol, Dls, Dos, tEin_scale, source_plane_FOV, arc_angle, kappa_micro, ml_position, ml_mass)


#### Writing the relevant input parametrs
plots.parameters(ml_position, ml_mass, tEin_scale, kappa_micro)


#### End of the program
