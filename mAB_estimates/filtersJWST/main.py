import numpy as np
import matplotlib.pyplot as plt


fn1 = 'F444W'
fn2 = '_NRC_and_OTE_ModAB_mean'
fn3 = '.txt'

filename = fn1 + fn2 + fn3
wavelength, throughput = np.loadtxt(filename, unpack=True, comments='#')

norm = 0.
for i in range(len(wavelength)-1):
    dx = wavelength[i+1] - wavelength[i]
    fx = 0.5*(throughput[i+1] + throughput[i])
    norm = norm + fx*dx

norm_throughput = throughput/norm

fname = 'normalized_' + fn1 + fn2 + fn3
np.savetxt(fname, list(zip(wavelength/1.e6, 299792458*1.e6/wavelength, norm_throughput)), delimiter='\t')

norm = 0.
for i in range(len(wavelength)-1):
    dx = wavelength[i+1] - wavelength[i]
    fx = norm_throughput[i]
    norm = norm + fx*dx

plt.plot(wavelength, norm_throughput)

plt.show()
