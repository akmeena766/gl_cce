import numpy as np
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


import input
from constants import *

import os
import time

# Setting up the home directory for to save plots and data
home_dir = os.getcwd()

# Path to the folder
path = home_dir + '/output/'



def plot_critical(tEin_scale, theta_x, theta_y, mu_lens, ml_position, ml_mass, theta_boundary):
    np.savetxt(path + 'microlens.dat',np.column_stack((ml_position, ml_mass)), fmt='%8.8f', delimiter='\t')
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    fig.subplots_adjust(left=.14, bottom=.1, right=.97, top=.9)
    
    cmap = 'cubehelix'  # the color map to use for the plot    
    plot_units = r'Angular position (in units of $\theta_{\rm{S,Ein}}$)'

    ax.set_xlabel(plot_units)
    ax.set_ylabel(plot_units)
    ax.set_aspect('equal')

    xmin, xmax = np.min(theta_x),np.max(theta_x)
    ymin, ymax = np.min(theta_y),np.max(theta_y)
    
    cs1 = ax.imshow(abs(mu_lens), origin = 'lower', extent = (xmin, xmax, ymin, ymax), cmap=cmap, norm=mpl.colors.LogNorm(vmin=1e0, vmax=1e3), aspect='auto')
    
    # Plotting the color bar
    cbar = plt.colorbar(cs1, pad=0.01)
    cbar.ax.set_ylabel('Magnification ($\mu$)', rotation=270, labelpad=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    # Plotting the point mass lenses
    ax.scatter(ml_position[:,0], ml_position[:,1], s = ml_mass, c='r', edgecolor='none')
    
    patch1 = patches.Rectangle((-theta_boundary[0], -theta_boundary[1]),
                                2*theta_boundary[0], 2*theta_boundary[1], linewidth=1, edgecolor='b', facecolor='none')

    ax.add_patch(patch1)
    
    ax.grid(alpha=0.5)

    # Adding twin x-axis in micro-arcseconds
    ax2 = ax.twiny()
    ax2.set_xticks( ax.get_xticks() )
    ax2.set_xbound(ax.get_xbound())
    ax2.set_xticklabels(['{:.0f}'.format(x*tEin_scale*ANGLE_ARCSEC/ANGLE_MICROARCSEC) for x in ax.get_xticks()])
    ax2.set_xlabel('Micro-arcseconds')

    plt.savefig(path + 'critcal_curve.jpeg', dpi=600)
    plt.close(fig)

    
def plot_caustics(tEin_scale, beta_x, beta_y, mu_source, source_plane_FOV):
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    fig.subplots_adjust(left=.13, bottom=.1, right=.97, top=.9)

    cmap = 'cubehelix'  # the color map to use for the plot

    plot_units = r'Angular position (in units of $\theta_{\rm{S,Ein}}$)'

    ax.set_xlabel(plot_units)
    ax.set_ylabel(plot_units)
    ax.set_aspect('equal')

    xmin, xmax = np.min(beta_x),np.max(beta_x)
    ymin, ymax = np.min(beta_y),np.max(beta_y)
    
    cs1 = ax.imshow(mu_source + 1, origin = 'lower', extent = (xmin, xmax, ymin, ymax), cmap=cmap, norm=mpl.colors.LogNorm(vmin=1e0, vmax=1e3), aspect='auto')
    
    # Plotting the color bar
    cbar = plt.colorbar(cs1, pad=0.01)
    cbar.ax.set_ylabel('Magnification ($\mu$)', rotation=270, labelpad=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    
    # Adding the relevant patch
    #patch = patches.Rectangle((-source_plane_FOV/2, -source_plane_FOV/2),
    #                            source_plane_FOV, source_plane_FOV, linewidth=1, edgecolor='r', facecolor='none')
    #ax.add_patch(patch)
    
    ax.grid(alpha=0.5)
    
    # Adding twin x-axis in micro-arcseconds
    ax2 = ax.twiny()
    ax2.set_xticks( ax.get_xticks() )
    ax2.set_xbound(ax.get_xbound())
    ax2.set_xticklabels(['{:.2f}'.format(x*tEin_scale*ANGLE_ARCSEC/ANGLE_MICROARCSEC) for x in ax.get_xticks()])
    ax2.set_xlabel('Micro-arcseconds')
    
    if input.method == 'IRS':                        
        plt.savefig(path + 'caustics_IRS.jpeg', dpi=600)
    
    plt.close(fig)


def plot_light_curve(time_steps, light_curve, size, shear_parallel):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    fig.subplots_adjust(left=.1, bottom=.11, right=.94, top=.985)

    plot_x_units = 'Time (in years)'
    plot_y_units = r'${\rm Log}_{10}(\mu)$'

    ax.set_xlabel(plot_x_units)
    ax.set_ylabel(plot_y_units)
    
    ax.set_xlim(0, input.light_curve_yr)
    
    ax.plot(time_steps, np.log10(light_curve + 1), c='k')
    ax.scatter(time_steps, np.log10(light_curve + 1), c='b')
    
    ax.grid(alpha=0.5)
    

    if input.method == 'IRS':
        size = '{:1.0e}'.format(size)
        if shear_parallel == True:
            plt.savefig(path + 'light_curve_parallel_IRS_' + size + '.jpeg', dpi=600)
            np.savetxt(path + 'light_curve_parallel_IRS_' + size + '.dat',list(zip(time_steps,light_curve)))
        else: 
            plt.savefig(path + 'light_curve_vertical_IRS_' + size + '.jpeg', dpi=600)
            np.savetxt(path + 'light_curve_vertical_IRS_' + size + '.dat',list(zip(time_steps,light_curve)))
    
    if input.method == 'ABM':
        size = '{:1.0e}'.format(input.source_size)
        if shear_parallel == True:
            plt.savefig(path + 'light_curve_parallel_ABM_' + size + '.jpeg', dpi=600)
            np.savetxt(path + 'light_curve_parallel_ABM_' + size + '.dat',list(zip(time_steps,light_curve)))
        else: 
            plt.savefig(path + 'light_curve_vertical_ABM_' + size + '.jpeg', dpi=600)
            np.savetxt(path + 'light_curve_vertical_ABM_' + size + '.dat',list(zip(time_steps,light_curve)))
    
    plt.close(fig)


def parameters(ml_position, ml_mass, tEin_scale, stelDen):
    # Saving the microlenses
    np.savetxt(path + 'microlens.dat',np.column_stack((ml_position*tEin_scale, ml_mass)), fmt='%8.8f', delimiter='\t')

    
    # Saving the input parameters
    file = open(path + 'parameters.txt', 'w')

    file.write('#### Cosmology ####\n')
    file.write('OM_0 = %3.4f\n' % input.OM_0)
    file.write('OB_0 = %3.4f\n' % input.OB_0)
    file.write('OL_0 = %3.4f\n' % input.OL_0)
    file.write('\n')
    file.write('H0 = %3.4f\n' % input.H0)

    file.write('\n')
    file.write('\n')

    file.write('#### Source ####\n')
    file.write('source_size = %3.3e\n' % input.source_size)
    file.write('vt = %3.4f\n' % input.vt)
    file.write('light_curve_yr = %3.4f\n' % input.light_curve_yr)

    file.write('\n')
    file.write('\n')

    file.write('#### Macrolens ####\n')
    file.write('kappa_macro = %3.4f\n' % input.kappa_macro)
    file.write('gamma_macro = %3.4f\n' % input.gamma_macro)
    file.write('kappa_frac = %3.4f\n' % input.kappa_frac)
    file.write('Stellar density = %5.4f (in MASS_SUN/PC^2)\n' % stelDen)

    file.write('\n')
    file.write('\n')

    file.write('#### Microlens population ####\n')
    file.write('lower_limit = %3.4f (in MASS_SUN)\n' % input.lower_limit)
    file.write('upper_limit = %3.4f (in MASS_SUN)\n' % input.upper_limit)
    file.write('alpha = %3.4f\n' % input.alpha)
    file.write('M_scale = %3.4f\n' % input.M_scale)

    file.write('\n')
    file.write('\n')

    file.write('#### Microlens simulation ####\n')
    file.write('user_seed = %3.4f\n' % input.user_seed)

    file.write('\n')

    file.write('method = %s\n' % input.method)
    if input.method == 'IRS':
        file.write('beta_res = %3.4f\n' % input.beta_res)
        file.write('rays_per_pixel = %3.4f\n' % input.rays_per_pixel)

    if input.method == 'ABM':
        file.write('n_pixels = %3.4f\n' % input.n_pixels)
        file.write('eta_ratio = %3.4f\n' % input.eta_ratio)
        file.write('n_steps = %3.4f\n' % input.n_steps)
        file.write('beta_0 = %3.4f\n' % input.beta_0)
        file.write('kernel_flag = %s\n' % input.kernel_flag)

