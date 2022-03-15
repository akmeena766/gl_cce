from astropy.cosmology import LambdaCDM
import astropy.units as units

import input as input

OM_0 = input.OM_0
OB_0 = input.OB_0
OL_0 = input.OL_0

H0 = input.H0

cosmo = LambdaCDM(H0=H0, Om0 = OM_0, Ob0 = OB_0, Ode0 = OL_0)

# Calculating Comoving Distance in meters
def comDist(z):
    return((cosmo.comoving_distance(z).to(units.m)).value)

# Calculating Angular Diameter Distance in meters
def angDist(z1,z2):
    D = cosmo.angular_diameter_distance_z1z2(z1,z2).to(units.m) 
    D = D.value             
    return(D)
