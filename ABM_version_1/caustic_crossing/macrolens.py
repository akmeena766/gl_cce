import numpy as np

import input
from constants import *
import cosmology as cosmo
import helper_functions as hf

rc = 50.*DIST_KPC		# Core radius in physical units
eps = 0.1			# Ellipticity of the lens
vdisp = 1000*1.e3 		# velocity dispersion of the lens


def macrolens_kgamma_NSIE(Dol, Dls, Dos, angle, kappa_micro, xx, yy):
    # Rotating the image to align it along the y-axis
    x1 = xx*np.cos(angle) - yy*np.sin(angle) 
    x2 = xx*np.sin(angle) + yy*np.cos(angle)
    
    # Potential multipication factor     
    psi_0 = 4.*np.pi*vdisp*vdisp*Dls/CONST_C/CONST_C/Dos
    
    # Core radius (in arcsec)
    tc = rc/Dol
    
    # Center of the lens plane patch (on the macro critical line)
    xc, yc = 0.*ANGLE_ARCSEC, -6.486965218141406*ANGLE_ARCSEC
    
    
    t1, t2 = x1*ANGLE_ARCSEC + xc, x2*ANGLE_ARCSEC + yc

    dfac = tc**2 + (1-eps)*t1**2 + (1+eps)*t2**2
    
    # Jacobian components
    psi_11 = psi_0*((1-eps)/np.sqrt(dfac) - (1-eps)**2*t1**2/dfac**(3/2))
    psi_22 = psi_0*((1+eps)/np.sqrt(dfac) - (1+eps)**2*t2**2/dfac**(3/2))
    psi_12 = -psi_0*((1-eps)*(1+eps)*t1*t2/dfac**(3/2))
    
    kappa = (psi_11 + psi_22)/2 #- kappa_micro
    gamma_1 = (psi_11 - psi_22)/2
    gamma_2 = psi_12
    
    gamma = np.sqrt(gamma_1**2 + gamma_2**2)
    
    return(kappa, gamma)
