import numpy as np
import random
from random import choices

import multiprocessing
from joblib import Parallel, delayed

import matplotlib as mpl
import matplotlib.patches
import matplotlib.colors
import matplotlib.pyplot as plt

from constants import *
import input as input
import macrolens as sl

#### Critical density
def critical_density(Dol, Dls, Dos):
    return(CONST_C**2.*Dos*Dol**2.*ANGLE_ARCSEC**2./(4.*np.pi*CONST_G*Dol*Dls*MASS_SUN))


#### kappa_micro to Stellar density
def kappaMicro_to_stelDen(Dol, critDen):
    return(input.kappa_micro*critDen*DIST_PC**2/(Dol**2*ANGLE_ARCSEC**2.))


#### Stellar density to kappa_micro
def stelDen_to_kappaMicro(stelDen, Dol, critDen):
    return(stelDen*(Dol**2*ANGLE_ARCSEC**2.)/critDen/DIST_PC**2)


#### Einstein angle for point mass lens (input mass in MASS_SUN)
def Einstein_Angle(Dol, Dls, Dos, M):
    return(np.sqrt(4.*CONST_G*M*MASS_SUN*Dls/(CONST_C**2.*Dol*Dos))/ANGLE_ARCSEC)


#### Source plane squre region size (in arcseconds)
def source_plane_FOV(Dos):
    return(input.vt*input.light_curve_yr*year2sec/Dos/ANGLE_ARCSEC)
    

#### Image plane squre region size (in arcseconds)
def image_plane_FOV(source_plane_FOV):
    return([input.epsilon_ml_draw*np.sqrt(input.mu_t**2. + input.mu_r**2.)*source_plane_FOV/2.,
            input.epsilon_ml_draw*np.sqrt(input.mu_t**2. + input.mu_r**2.)*source_plane_FOV/2.])


#### Generating the microlens population
def generate_microlens_population(kappa_micro, image_plane_FOV, critDen, arc_angle):
    lower_limit = input.lower_limit
    upper_limit = input.upper_limit
    alpha = input.alpha
    
    kappa_micro = abs(kappa_micro)		  			# Microlens convergence
    
    stelden = kappa_micro*critDen					# Stellar Density (in MASS_SUN/arcsec^2)

    image_plane_area = 4*image_plane_FOV[0]*image_plane_FOV[1]		# Image plane area to put microlenses (arcsec^2)
    
    total_mass = kappa_micro*critDen*image_plane_area			# Total mass of the microlenses (in MASS_SUN)

    max_draw = int(total_mass/lower_limit + 0.5)				# Maximum possible number of microlenses
    
    microlens_mass = np.arange(lower_limit, upper_limit, 1e-4)
    imf = microlens_mass**(-alpha)					# IMF function: \propto m**(-alpha)
    
    # Generating microlens masses
    itr = 1
    total_drawn_mass = 0 
    ml_mass = np.array([])
    while total_drawn_mass < total_mass:
          ndraw = int(max_draw/10**itr)

          if ndraw < 100:
             for i in range(100):
                 draw_mass = choices(microlens_mass, imf, k=1)
                 total_drawn_mass = total_drawn_mass + sum(draw_mass)
                 if(total_drawn_mass < total_mass):   
                     ml_mass = np.append(ml_mass,draw_mass[0])
                 else:
                     break
          else:
              draw_mass = choices(microlens_mass, imf, k=ndraw)
              total_drawn_mass = total_drawn_mass + sum(draw_mass)

              if(total_drawn_mass < total_mass):   
                  ml_mass = np.append(ml_mass,draw_mass)
              else:
                  itr = itr + 1
                  total_drawn_mass = total_drawn_mass - sum(draw_mass)
              
    # Generarting microlens postions    
    position_x = image_plane_FOV[0]*(2*np.random.rand(len(ml_mass)).reshape(len(ml_mass), 1) - 1)
    position_y = image_plane_FOV[1]*(2*np.random.rand(len(ml_mass)).reshape(len(ml_mass), 1) - 1)

    ml_position = np.hstack((position_x, position_y))
    return(ml_position, ml_mass)


#### parallel function to calculate magnification
def magnification_function(i,theta_x, theta_y_value, tEin_scale, Dol, Dls, Dos, arc_angle, kappa_micro, ml_x, ml_y, ml_mass):
    theta_y = np.full((theta_x.shape[0],1), theta_y_value)

    dx = theta_x - ml_x
    dy = theta_y - ml_y
    
    psi_11 = ((dy**2 - dx**2)/(dx**2 + dy**2)**2) @ (ml_mass/input.M_scale)  	# size(l,1) via matrix multipication
    psi_22 = -psi_11  					    			# size(l,1) via matrix multipication
    psi_12 = (-2*dx*dy/(dx**2 + dy**2)**2) @ (ml_mass/input.M_scale) 	        # size(l,1) via matrix multipication

    gamma1 = (psi_11 - psi_22)/2
    gamma2 = psi_12
    
    ks, gs = sl.macrolens_kgamma_NSIE(Dol, Dls, Dos, arc_angle, kappa_micro, theta_x*tEin_scale, theta_y*tEin_scale)
    
    mag = (1. - ks)**2 - (gs + gamma1)**2 - gamma2**2
    return(1/mag)


#### Calculating the image plane grid magnification values (only for IRS)
def critical_curves(tEin_scale, Dol, Dls, Dos, theta_boundary, arc_angle, kappa_micro, ml_position, ml_mass):
    #### Pixelizing the lens plane (By default, always grid a square region)
    square = True

    #### If the given resolution is very small
    if input.beta_res < 1000:
        resolution = 1000
    else:
        resolution = input.beta_res
        
    if square == True:
        theta_x = np.arange(-theta_boundary[1], theta_boundary[1], 2.*theta_boundary[1]/resolution)
        theta_y = theta_x
    else:
        theta_x = np.arange(-theta_boundary[0], theta_boundary[0], 2.*theta_boundary[0]/resolution)
        theta_y = np.arange(-theta_boundary[1], theta_boundary[1], 2.*theta_boundary[1]/resolution)

    diff_len = len(theta_x) - len(theta_y)

    if diff_len > 0: theta_x = theta_x[:-diff_len]
    if diff_len < 0: theta_y = theta_y[:+diff_len]
    
    theta_x = theta_x.reshape(theta_x.shape[0],1)   # size(l,1)
    theta_y = theta_y.reshape(theta_y.shape[0],1)   # size(l,1)
    
    ml_mass = ml_mass.reshape((ml_mass.shape[0], 1))  # size(n,1)
    
    ml_x = ml_position[:, 0].reshape((1, ml_mass.shape[0]))  # size(1,n)
    ml_y = ml_position[:, 1].reshape((1, ml_mass.shape[0]))  # size(1,n)
    
    num_cores = multiprocessing.cpu_count()
       
    def_tmp = Parallel(n_jobs=num_cores)\
    (delayed(magnification_function)(i, theta_x, theta_y[i,0], tEin_scale, Dol, Dls, Dos, 
                                     arc_angle, kappa_micro, ml_x, ml_y, ml_mass) for i in range(theta_x.shape[0]))

    def_tmp = np.concatenate(def_tmp, axis=0)
    
    magnification = np.reshape(def_tmp, (theta_x.shape[0],theta_y.shape[0]))

    return(theta_x, theta_y, magnification)
    
    
#### Tesellating the image plane (relevant for ABM and IPM)
def tessellate_area(theta_boundary, offset, num_of_tiles, ret_vertices=True, ret_centers=True):
    # x and y axis half limits
    x_lim, y_lim = theta_boundary[0], theta_boundary[1]
    
    # Number of tiles per row
    num_of_tiles_per_row = int(np.sqrt(num_of_tiles))
    
    # Retrieveing vertices
    if ret_vertices:
        x = np.linspace(-x_lim, x_lim, num_of_tiles_per_row + 1, endpoint=True) + offset[0]
        y = np.linspace(-y_lim, y_lim, num_of_tiles_per_row + 1, endpoint=True) + offset[1]

        grid_x, grid_y = np.meshgrid(x, y)
        
        true_num_tiles = grid_x.shape[0] * grid_x.shape[1]

        vertices = np.vstack((grid_x.reshape((1, true_num_tiles)), grid_y.reshape((1, true_num_tiles)))).T

    # Retrieveing center values
    if ret_centers:
        step_x = x_lim * 2 / (num_of_tiles_per_row + 1)
        step_y = y_lim * 2 / (num_of_tiles_per_row + 1)

        x = np.linspace(-x_lim + step_x / 2, x_lim - step_x / 2, num_of_tiles_per_row, endpoint=True) + offset[0]
        y = np.linspace(-y_lim + step_y / 2, y_lim - step_y / 2, num_of_tiles_per_row, endpoint=True) + offset[1]
        
        grid_x, grid_y = np.meshgrid(x, y)
        
        true_num_tiles = grid_x.shape[0] * grid_y.shape[1]
        
        centers = np.vstack((grid_x.reshape((1, true_num_tiles)), grid_y.reshape((1, true_num_tiles)))).T
        
        #### Retrieveing vertices and center values
        if ret_vertices:
            return vertices, centers
        else:
            return centers
    else:
        return vertices


#### Shooting rays from image to source plane
def img2src(tEin_scale, Dol, Dls, Dos, theta, arc_angle, kappa_micro, ml_position, ml_mass):
    n = ml_mass.shape[0]
    l = theta.shape[0]
    
    theta_x = theta[:, 0].reshape((l, 1))  # size(l,1)
    theta_y = theta[:, 1].reshape((l, 1))  # size(l,1) 
    
    ml_x = ml_position[:, 0].reshape((1, n))  # size(1,n)
    ml_y = ml_position[:, 1].reshape((1, n))  # size(1,n)
    ml_mass = ml_mass.reshape((n, 1))
    
    dx = theta_x - ml_x
    dy = theta_y - ml_y
    
    alpha_x = (dx/(dx**2 + dy**2 + 1.e-15)) @ (ml_mass/input.M_scale)   # size(l,1) via matrix multipication
    alpha_y = (dy/(dx**2 + dy**2 + 1.e-15)) @ (ml_mass/input.M_scale)   # size(l,1) via matrix multipication
    
    alpha = np.hstack((alpha_x.reshape((alpha_x.shape[0], 1)), 
                       alpha_y.reshape((alpha_y.shape[0], 1))))  # size(l,2)
    
    ks, gs = sl.macrolens_kgamma_NSIE(Dol, Dls, Dos, arc_angle, kappa_micro, theta_x*tEin_scale, theta_y*tEin_scale)
    
    beta_x = theta_x*(1- ks - gs) - alpha_x
    beta_y = theta_y*(1- ks + gs) - alpha_y
    
    beta = np.hstack((beta_x, beta_y))
                           
    return(beta)
    
    
#### Magnification binning in the source plane (relevant for IRS)
def mag_binning(beta, weight, boundary, beta_res = input.beta_res, beta_offset=[0,0]):
    min_x, max_x = -boundary[0]/2 + beta_offset[0], boundary[0]/2 + beta_offset[0]
    min_y, max_y = -boundary[1]/2 + beta_offset[1], boundary[1]/2 + beta_offset[1]
    
    mu_tot = np.zeros((beta_res, beta_res))

    if 4*beta_res**2 > 8*10**9:  # if file size is larger than 4GB
        return 0

    grid_vec_x, dx = np.linspace(min_x, max_x, beta_res, endpoint=False, retstep=True)
    grid_vec_y, dy = np.linspace(min_y, max_y, beta_res, endpoint=False, retstep=True)

    beta_grid_x, beta_grid_y = np.meshgrid(grid_vec_x, grid_vec_y)

    if beta.shape[0] > 0:  # if beta list is not empty
        beta_x = beta[:, 0].reshape((beta.shape[0], 1))  # size(l,1)
        beta_y = beta[:, 1].reshape((beta.shape[0], 1))  # size(l,1)
        
        mu_x_idx = np.squeeze(np.floor_divide(beta_x - min_x, dx).astype(int))
        mu_y_idx = np.squeeze(np.floor_divide(beta_y - min_y, dy).astype(int))
        
        for i in range(beta.shape[0]):
            if beta_res > mu_x_idx[i] >= 0 <= mu_y_idx[i] < beta_res:
                mu_tot[mu_x_idx[i], mu_y_idx[i]] += weight
    return beta_grid_x, beta_grid_y, mu_tot.T/(dx*dy)


