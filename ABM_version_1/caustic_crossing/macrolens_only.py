import numpy as np
import matplotlib.pyplot as plt

import input
from constants import *
import cosmology as cosmo
import helper_functions as hf

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig.subplots_adjust(left=.1, bottom=.1, right=.97, top=.9)

# Observer, Lens and Source redshifts
zO = 0.0
zL = 0.5
zS = 1.5


# Relevant Angular Diameter Distances (ADD)
Dol = cosmo.angDist(zO, zL)
Dls = cosmo.angDist(zL, zS)
Dos = cosmo.angDist(zO, zS)

# Critical density (in MASS_SUN/arcsec^2)
critDen = hf.critical_density(Dol,Dls,Dos)

rc = 50.*DIST_KPC		# Core radius (in meters)
eps = 0.1			# Ellipticity of the lens
vdisp = 1000*1.e3 		# velocity dispersion of the lens

# Factor in multipication to the potential
psi_0 = 4.*np.pi*vdisp*vdisp*Dls/CONST_C/CONST_C/Dos

tc = rc/Dol			# Core radius (in arcsec)

def nsie_kgamma(t1, t2):
    dfac = tc**2 + (1-eps)*t1**2 + (1+eps)*t2**2
    
    # potential function
    psi = psi_0*np.sqrt(dfac)
    
    # Deflection vector components
    psi_1 = psi_0*(t1*(1-eps)/np.sqrt(dfac))
    psi_2 = psi_0*(t2*(1+eps)/np.sqrt(dfac))
    
    # Jacobian components
    psi_11 = psi_0*((1-eps)/np.sqrt(dfac) - (1-eps)**2*t1**2/dfac**(3/2))
    psi_22 = psi_0*((1+eps)/np.sqrt(dfac) - (1+eps)**2*t2**2/dfac**(3/2))
    psi_12 = -psi_0*((1-eps)*(1+eps)*t1*t2/dfac**(3/2))
    
    kappa = (psi_11 + psi_22)/2
    gamma_1 = (psi_11 - psi_22)/2
    gamma_2 = psi_12
    
    gamma = np.sqrt(gamma_1**2 + gamma_2**2)
    
    detA = ((1-kappa)**2 - gamma**2)
    return(kappa, gamma, detA, gamma_1, gamma_2)
    
def nsie_alpha(t1, t2):
    dfac = tc**2 + (1-eps)*t1**2 + (1+eps)*t2**2
    
    # potential function
    psi = psi_0*np.sqrt(dfac)
    
    # Deflection vector components
    psi_1 = psi_0*(t1*(1-eps)/np.sqrt(dfac))
    psi_2 = psi_0*(t2*(1+eps)/np.sqrt(dfac))
    
    b1 = t1 - psi_1
    b2 = t2 - psi_2
    return(b1, b2)

xmin, xmax = -20, 20
ymin, ymax = -20, 20

x1, x2 = np.meshgrid(np.arange(xmin*ANGLE_ARCSEC, xmax*ANGLE_ARCSEC, (xmax-xmin)*ANGLE_ARCSEC/1000), 
                     np.arange(ymin*ANGLE_ARCSEC, ymax*ANGLE_ARCSEC, (ymax-ymin)*ANGLE_ARCSEC/1000))

kval, gval, detA, g1, g2 = nsie_kgamma(x1, x2)


c1 = ax2.contour(x1/ANGLE_ARCSEC,x2/ANGLE_ARCSEC,detA, [0], colors='red', zorder = 1)
crit_xt, crit_yt, dis = np.array([]), np.array([]), np.array([])
for i in range(len(c1.allsegs[0])):
    for j in range(len(c1.allsegs[0][i])):
        crit_xt = np.append(crit_xt, c1.allsegs[0][i][j,0])
        crit_yt = np.append(crit_yt, c1.allsegs[0][i][j,1])

xx, yy = 0, -6
dis = np.sqrt((crit_xt - xx)**2 + (crit_yt - yy)**2)
index = np.where(dis == dis.min())

xf, yf = crit_xt[index[0][0]], crit_yt[index[0][0]]
print(xf, yf)

kf, gf, af, g1f, g2f = nsie_kgamma(xf*ANGLE_ARCSEC, yf*ANGLE_ARCSEC)
print(g1f, g2f)


for i in range(len(crit_xt)):
    dis = np.sqrt((crit_xt[i] - xf)**2 + (crit_yt[i] - yf)**2)
    if dis < 1: 
        ax2.scatter(crit_xt[i], crit_yt[i], c='g', s=2, zorder = 2)
        beta_x, beta_y = nsie_alpha(crit_xt[i]*ANGLE_ARCSEC, crit_yt[i]*ANGLE_ARCSEC)
        ax1.scatter(beta_x/ANGLE_ARCSEC, beta_y/ANGLE_ARCSEC, c='g', s=2, zorder = 2)

ax2.set_xlabel('arcsec')
ax2.set_ylabel('arcsec')

ax2.grid(alpha = 0.5)



# Source plane plot
beta_x, beta_y = nsie_alpha(crit_xt*ANGLE_ARCSEC, crit_yt*ANGLE_ARCSEC)
ax1.scatter(beta_x/ANGLE_ARCSEC, beta_y/ANGLE_ARCSEC, c='r', s=1, zorder = 1)

ax1.set_xlim(-10,10)
ax1.set_ylim(-10,10)
ax1.set_xlabel('arcsec')
ax1.set_ylabel('arcsec')

ax1.grid(alpha = 0.5)
plt.show()
