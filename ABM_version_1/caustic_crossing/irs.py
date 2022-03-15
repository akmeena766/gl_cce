import numpy as np
import time

import multiprocessing
from joblib import Parallel, delayed

import input
import plots
import helper_functions as hf
import light_curve as lc


def irs_main(Dol, Dls, Dos, tEin_scale, source_plane_FOV, arc_angle, kappa_micro, ml_position, ml_mass):
    #### Region in the image plane to be grided (in arcsec)
    theta_boundary = [input.xbox - input.buff, input.ybox - input.buff]  
                      

    #### Region in image plane to be grided (in units of tEin_scale)
    theta_boundary = theta_boundary/tEin_scale

    #### Region in source plane to be grided (in units of tEin_scale)
    source_plane_FOV = [source_plane_FOV, source_plane_FOV]/tEin_scale + [0, 0]
    
    
    #### Critical curve generation in the lens plane
    if input.plot_critical:
        theta_x, theta_y, mu_lens = hf.critical_curves(tEin_scale, Dol, Dls, Dos, theta_boundary, arc_angle, kappa_micro, ml_position, ml_mass)
        plots.plot_critical(tEin_scale, theta_x, theta_y, mu_lens, ml_position, ml_mass, theta_boundary)
    
    
    #### Caustic curve generation in the source plane
    if input.plot_caustics == True:
        beta_x, beta_y, mu_source = irs_caustics_random_shooting(tEin_scale, Dol, Dls, Dos, theta_boundary, 
                                                                 source_plane_FOV, arc_angle, kappa_micro, ml_position, ml_mass)
        plots.plot_caustics(tEin_scale, beta_x, beta_y, mu_source, source_plane_FOV)
    
    
        #### Source size values for which we will construct the light curves
        source_size = [1.e11, 1.e12, 1.e13]
    
        #### Constructing the light curves
        for i in range(len(source_size)):
            #### Generating parallel light curve
            if input.plot_light_curve_parallel == True:
               time_steps, light_curve_values = lc.light_curve(Dos, source_plane_FOV[1], source_size[i], tEin_scale, mu_source, shear_parallel = True)
               plots.plot_light_curve(time_steps, light_curve_values, source_size[i], shear_parallel = True)
    
            #### Generating vertical light curve
            if input.plot_light_curve_vertical == True:
                time_steps, light_curve_values = lc.light_curve(Dos, source_plane_FOV[1], source_size[i], tEin_scale, mu_source, shear_parallel = False)
                plots.plot_light_curve(time_steps, light_curve_values, source_size[i], shear_parallel = False)
    


def irs_caustics_random_shooting(tEin_scale, Dol, Dls, Dos, theta_boundary, source_plane_FOV, arc_angle, kappa_micro, ml_position, ml_mass):
    #### Total number of pixels in the source plane
    num_of_pixel = input.beta_res**2.
    
    #### Total number of rays to shoot from image to source plane
    num_of_rays = num_of_pixel*input.rays_per_pixel
    
    #### The maximum number of rays to shoot at a time (to avoid memory issues)
    l_tmp = int(input.max_memory*1.e9/ml_mass.shape[0]/8.)

    #### Total number of runs
    n_runs = max(int(num_of_rays/l_tmp),1)
    
    #### Area associated with each ray in lens plane (in units of tEin_scale^2)
    area_per_ray = 4.*theta_boundary[0]*theta_boundary[1]/num_of_rays
    
    #### Initialize the source plane grid
    mu_grid = np.zeros((input.beta_res, input.beta_res))
    
    num_cores = multiprocessing.cpu_count()
    start_time = time.time()    

    mu_grid_temp_array = Parallel(n_jobs=num_cores, require='sharedmem')\
        (delayed(parallel_irs)(i, num_of_rays, n_runs, start_time, mu_grid, l_tmp, area_per_ray, source_plane_FOV, 
         tEin_scale, Dol, Dls, Dos, theta_boundary, arc_angle, kappa_micro, ml_position, ml_mass) for i in range(n_runs))


    if n_runs * l_tmp < num_of_rays:  # if some values are left
        # Drawing images locations
        theta = random_image_draw(int(num_of_rays-(n_runs*l_tmp)), theta_boundary[0], theta_boundary[1])
        
        # Calculating locations of sources and corresponding magnitudes
        beta = hf.img2src(tEin_scale, Dol, Dls, Dos, theta, arc_angle, kappa_micro, ml_position, ml_mass)
        
        # Binning sources magnification
        beta_grid_h, beta_grid_v, mu_grid_temp = hf.mag_binning(beta, area_per_ray, source_plane_FOV)
        mu_grid += mu_grid_temp
    else:
        beta = np.ones((2, 2))  # Just so that the next line can run smoothly and return beta_grid_h and beta_grid_v
        beta_grid_h, beta_grid_v, mu_grid_temp = hf.mag_binning(beta, area_per_ray, source_plane_FOV)

    return beta_grid_h, beta_grid_v, mu_grid



#### Parellel part of the IRS
def parallel_irs(i, num_of_rays, n_runs, start_time, mu_grid, l_tmp, area_per_ray, source_plane_FOV, 
                 tEin_scale, Dol, Dls, Dos, theta_boundary, arc_angle, kappa_micro, ml_position, ml_mass):
    # Drawing images locations theta
    theta = random_image_draw(l_tmp, theta_boundary[0], theta_boundary[1])
    
    # Calculating locations of sources and corresponding magnitudes
    beta = hf.img2src(tEin_scale, Dol, Dls, Dos, theta, arc_angle, kappa_micro, ml_position, ml_mass)
    
    # Binning sources magnification
    beta_grid_x, beta_grid_y, mu_grid_temp = hf.mag_binning(beta, area_per_ray, source_plane_FOV)
    mu_grid += mu_grid_temp

    temp_t = time.time() - start_time
    if i % (max(1, int(n_runs / 4000))) == 0:
        print('Inverse Ray Shooting: Finished ' + str(round((i + 1) * 100 / n_runs, 3)) + '% in ' + str(round(temp_t)) +
              's; ~' + str(round(temp_t * (n_runs / (i + 1) - 1))) + 's remaining')
    return(0)


#### Random position in the image plane to shoot
def random_image_draw(num_of_img, h_lim, v_lim):
    theta = np.hstack((h_lim * (2 * np.random.rand(num_of_img).reshape(num_of_img, 1) - 1),
                       v_lim * (2 * np.random.rand(num_of_img).reshape(num_of_img, 1) - 1)))
    return(theta)
    
