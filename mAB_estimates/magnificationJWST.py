import numpy as np
from astropy.cosmology import LambdaCDM
import astropy.units as units
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d


# Relevant constants
OM_0 = 0.3				# Present matter density
OB_0 = 0.0				# Present baryon density
OL_0 = 0.7				# Cosmological constant 

H0 = 70.				# Hubble constant

c = 29979245800.0			# Speed of light
h = 6.62607004e-27			# Planck constant
k = 1.38064852e-16			# Boltzmann constant
RADIUS_SUN = 695700.e5			# Radius of the Sun

DIST_PC = 3.0856775714409184e+18	# One parsec (in cm)

cosmo = LambdaCDM(H0=H0, Om0 = OM_0, Ob0 = OB_0, Ode0 = OL_0)


# luminosity distance
def lumDist(z):
    return(cosmo.luminosity_distance(z).to(units.cm).value)


# Calling the relevant filter
fn1 = './filtersJWST/normalized_'
fn2 = 'F200W'
fn3 = '_NRC_and_OTE_ModAB_mean'
fn4 = '.txt'
print(fn2)
filename = fn1 + fn2 + fn3 + fn4
wavelength, frequency, norm_throughput = np.loadtxt(filename, unpack=True, comments='#')


# Zero point magnitude for any band in AB system
zero_point_mag = 48.60


# Normalization constant
norm = 0
for i in range(len(frequency)-1):
    df = abs(frequency[i+1] - frequency[i])
    f = frequency[i]
    Qf = norm_throughput[i+1]
    norm = norm + Qf*df/f


aa, bb = [], []
# Absolute mangitude calculation
def abs_mag(data, zS, T, rS):
    source_flux = 0
    B = interp1d(data[:,0], data[:,1])
    for i in range(len(frequency)-1):
        df = abs(frequency[i+1] - frequency[i])
        f = frequency[i]
        Qf = norm_throughput[i]
        source_flux = source_flux + B(f)*Qf*df/f

    source_flux = source_flux*(1+zS)*np.pi*rS**2/(10.*DIST_PC)**2
    MAB = - 2.5*np.log10(source_flux/norm) - zero_point_mag
    return(MAB)



# Apparent magnitude calculation
def app_mag(data, frequency, zS, T, rS):
    source_flux = 0
    frequency = frequency*(1+zS)
    B = interp1d(data[:,0], data[:,1])
    for i in range(len(frequency)-1):
        df = abs(frequency[i+1] - frequency[i])
        f = frequency[i]
        Qf = norm_throughput[i]
        source_flux = source_flux + B(f)*Qf*df/f
        
    source_flux = source_flux*(1+zS)*np.pi*rS**2/(lumDist(zS))**2
    mAB = - 2.5*np.log10(source_flux/norm) - zero_point_mag
    return(mAB)


zS = 1.5

Temp = [5750, 30000, 45000, 4000, 12000]
rS = [1*RADIUS_SUN, 6.5*RADIUS_SUN, 12*RADIUS_SUN, 45*RADIUS_SUN, 79*RADIUS_SUN]

print('Temperature \t MAB \t\t mAB \t\t 27AB \t\t 29AB \t\t 31AB\n')
for i in range(len(Temp)):
    fname = './stellarSpectra/ckp00_'+ str(Temp[i]) +'.dat'
    data = np.loadtxt(fname, unpack=True)
    data = data.T

    M_AB = abs_mag(data, zS = 0., T = Temp[i], rS = rS[i])

    m_AB = app_mag(data, frequency, zS = zS, T = Temp[i], rS = rS[i])

    e_AB = 27
    mu1 = 10**((m_AB - e_AB)/2.5)  

    e_AB = 29
    mu2 = 10**((m_AB - e_AB)/2.5)  

    e_AB = 31
    mu3 = 10**((m_AB - e_AB)/2.5)  
    print('%f\t%f\t%f\t%4.1e\t\t%4.1e\t\t%4.1e'%(Temp[i], M_AB, m_AB, round(mu1,1), round(mu2,1), round(mu3,1)))
