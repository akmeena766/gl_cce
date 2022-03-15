import numpy as np
import time

import multiprocessing
from joblib import Parallel, delayed

import input
import plots
from constants import *
import helper_functions as hf


def abm_main(Dol, Dls, Dos, tEin_scale, source_plane_FOV, ml_position, ml_mass):
    # Region in the lens plane to be grided (in arcsec)
    theta_boundary = [input.epsilon_ml_shoot * input.mu_t * source_plane_FOV/2.,
                      input.epsilon_ml_shoot * input.mu_t * source_plane_FOV/2.]                              

    # Region in lens and source plane to be grided (in units of tEin_scale)
    theta_boundary = theta_boundary/tEin_scale
    source_plane_FOV = source_plane_FOV/tEin_scale
    
    # Initial source plane radial-boundary for the ABM algorithm
    beta_boundary_0 = input.beta_0*source_plane_FOV

    # Final source plane half boundary (in Einstein radii)
    delta_beta = input.source_size/Dos/ANGLE_ARCSEC/tEin_scale
    
    # eta: Factor by which the initial beta boundary will be refined in each iteration
    # #### around a source point
    eta = input.n_pixels*input.eta_ratio
    
    # Number of iterations (NOTE: it is a real number, not an integer)
    num_iter = 1 + np.log(beta_boundary_0/delta_beta)/np.log(eta)
    print(num_iter)
    # Source velocity (in Einstein_angle/second)
    vt_sEin = input.vt/Dos/ANGLE_ARCSEC/tEin_scale
    
    # Calculating the parallel light curve path
    if input.plot_light_curve_parallel:
        path_vec = np.linspace(-source_plane_FOV/2, source_plane_FOV/2, num = input.n_steps, endpoint=True)
        path_vec = np.vstack((path_vec, np.zeros((1,path_vec.shape[0])))).T
        

        # Calculating the parallel light curve
        light_curve = abm_light_curve(tEin_scale, Dol, Dls, Dos, path_vec, theta_boundary, beta_boundary_0, input.n_pixels, eta, 
                                  num_iter, delta_beta, ml_position, ml_mass, 'parallel')
         
        # Light curve time steps       
        time_steps = path2displacement(path_vec)*(1/vt_sEin/year2sec)
        
        # Plotting the light curve
        plots.plot_light_curve(time_steps, light_curve, 0, shear_parallel = True)


    # Calculating the vertical light curve
    if input.plot_light_curve_vertical:
        path_vec = np.linspace(-source_plane_FOV/2, source_plane_FOV/2, num = input.n_steps, endpoint=True)
        path_vec = np.vstack((np.zeros((1, path_vec.shape[0])), path_vec)).T
                
        # Calculating the vertical light curve
        light_curve = abm_light_curve(tEin_scale, Dol, Dls, Dos, path_vec, theta_boundary, beta_boundary_0, input.n_pixels, eta, 
                                  num_iter, delta_beta, ml_position, ml_mass, 'vertical')
        
        # Light curve time steps
        time_steps = path2displacement(path_vec)*(1/vt_sEin/year2sec)  
        
        # Plotting the light curve
        plots.plot_light_curve(time_steps, light_curve, 0, shear_parallel = False)



def abm_light_curve(tEin_scale, Dol, Dls, Dos, path_vec, theta_lim_0, beta_lim_0, n_pixels, eta, num_iter, delta_beta, ml_position, ml_mass, lcd):
    # Image plane area of each ray
    s_img = 4*theta_lim_0[0]*theta_lim_0[1]/n_pixels**(2*np.floor(num_iter) + 2)
    
    # Whether to apply Gaussian profile before creating light curve
    if input.kernel_flag:
        kernel = gaussian_kernel(3, 9)  # 9X9 matrix, corresponding to a range of +-3*sigma of a Gaussian kernel

    # Initializing empty light curve
    light_curve = []
    
    ### In the global iteration, we are getting relevant pixels for the whole light curve in one go.
    ### We are doing this for some iteration such that we do not run into RAM issues.
    if input.sub_light_curve_yr == 0:
        multi = input.n_steps
    else:
        multi = int(input.light_curve_yr/input.sub_light_curve_yr)
    
    for j in range(multi):
        start_time = time.time()
        step_size = int(input.n_steps/multi)
        path_vec1 = path_vec[j*step_size : (j+1)*step_size]
        if j == multi-1:
            path_vec1 = path_vec[j*step_size : ]

        new_beta_boundaries = beta_lim_0
        new_theta_boundaries = theta_lim_0
    
        all_thetas = np.zeros((1,2))
    
        for i in range(0, int(num_iter) + 1):
            if i < int(num_iter):
                new_beta_boundaries = beta_lim_0/eta**i
                new_theta_boundaries = [theta_lim_0[0]/n_pixels**i, theta_lim_0[1]/n_pixels**i]

                global_theta = split_and_identify(tEin_scale, Dol, Dls, Dos, new_theta_boundaries, all_thetas, new_beta_boundaries, 
                                                               path_vec1, n_pixels, ml_position, ml_mass, lcd, 'global')                                        

                all_thetas = global_theta
                
            else:
                new_theta_boundaries = [new_theta_boundaries[0]/n_pixels, new_theta_boundaries[1]/n_pixels]
            
                # Decrease the source plane boundary by a factor eta^n2  
                n1 = int(np.floor(num_iter))
                n2 = num_iter - n1
                new_beta_boundaries = new_beta_boundaries/eta**n2

                global_beta = split_and_identify(tEin_scale, Dol, Dls, Dos, new_theta_boundaries, all_thetas, new_beta_boundaries, 
                                                               path_vec1, n_pixels, ml_position, ml_mass, lcd, 'final')
    
        num_cores = multiprocessing.cpu_count()   
        # Binning the rays in the source plane
        tmp_light_curve = Parallel(n_jobs=num_cores, prefer='threads')\
                      (delayed(bin_rays)(i, path_vec1[i], global_beta, new_beta_boundaries, s_img, delta_beta) for i in range(len(path_vec1)))
        light_curve.append(tmp_light_curve)
        print('ABM Light Curve: ' + str(step_size*(j+1)) + ' points out of ' + str(input.n_steps) + ' are done in ' + 
        str(round(time.time()-start_time, 2)))
    light_curve = np.concatenate(light_curve)
    return(np.array(light_curve))    
    


def split_and_identify(tEin_scale, Dol, Dls, Dos, theta_lim, theta_offset, beta_lim, beta_offset, n_pixels, ml_position, ml_mass, lcd, itrv):
    num_cores = multiprocessing.cpu_count()  
    
    # Split the image plane region to nxn pixels
    theta_pixels = Parallel(n_jobs=num_cores)\
                   (delayed(hf.tessellate_area)(theta_lim, theta_offset[i], n_pixels**2, ret_vertices=False) for i in range(len(theta_offset)))
    theta_pixels = np.concatenate(theta_pixels, axis=0)

    # The maximum number of rays to shoot at a time 
    num_of_rays = len(theta_pixels)
    l_tmp = int(input.max_memory*1.e9/ml_mass.shape[0]/8.)

    # Total number of runs
    n_runs = max(int(num_of_rays/l_tmp),1)

    if n_runs > 1:
        tmp_beta_pixels = Parallel(n_jobs=num_cores)\
        (delayed(hf.img2src)(tEin_scale, Dol, Dls, Dos, theta_pixels[i*l_tmp : (i+1)*l_tmp], ml_position, ml_mass) for i in range(n_runs-1))
        
        tmp_beta_pixels.append(hf.img2src(tEin_scale, Dol, Dls, Dos, theta_pixels[(n_runs-1)*l_tmp : ], ml_position, ml_mass))
        
        beta_pixels = np.concatenate(tmp_beta_pixels)
    else:
        beta_pixels = hf.img2src(tEin_scale, Dol, Dls, Dos, theta_pixels, ml_position, ml_mass)    

    # Identify relevant pixels, that map to the source plane region of interest
    if input.point_based:
        rel_pixels = [np.nonzero(((beta_offset[step, 0] - beta_pixels[:, 0])**2 + (beta_offset[step, 1] - beta_pixels[:, 1])**2)
                      < beta_lim ** 2)[0] for step in range(len(beta_offset))]        
        rel_pixels = np.unique(np.concatenate(rel_pixels, axis=0))
        
        theta_pixels_out, beta_pixels_out = [], []
        for pixel_num in rel_pixels:
            theta_pixels_out.append(theta_pixels[pixel_num])
            beta_pixels_out.append(beta_pixels[pixel_num])  
        
        theta_pixels_out = np.reshape(theta_pixels_out, (len(theta_pixels_out), 2))  
        beta_pixels_out = np.reshape(beta_pixels_out, (len(theta_pixels_out), 2))  
        if itrv == 'final':
            return(beta_pixels_out)  
        else:
            return(theta_pixels_out) 
    else:
        beta_pixels = np.hstack((beta_pixels, theta_pixels))
        theta_pixels = np.zeros((1,2))
        
        if lcd == 'parallel':
            beta_pixels = beta_pixels[np.abs(beta_pixels[:, 1]) < beta_lim]
            rel_pixels = beta_pixels[(beta_pixels[:, 0] > beta_offset[0, 0] - beta_lim) & (beta_pixels[:, 0] < beta_offset[len(beta_offset)-1, 0] + beta_lim)]
        else:         
            beta_pixels = beta_pixels[np.abs(beta_pixels[:, 0]) < beta_lim]
            rel_pixels = beta_pixels[(beta_pixels[:, 1] > beta_offset[0, 1] - beta_lim) & (beta_pixels[:, 1] < beta_offset[len(beta_offset)-1, 1] + beta_lim)]
        
        beta_pixels = np.zeros((1,2))

        theta_pixels_out = np.vstack((rel_pixels[:, 2], rel_pixels[:, 3])).T
        beta_pixels_out  = np.vstack((rel_pixels[:, 0], rel_pixels[:, 1])).T
        if itrv == 'final':
            return(beta_pixels_out)  
        else:
            return(theta_pixels_out)



def gaussian_kernel(num_sigma, num_pixels):
    temp_x, temp_y = np.meshgrid(np.linspace(-num_sigma, num_sigma, num_pixels),
                                 np.linspace(-num_sigma, num_sigma, num_pixels))
    
    # Gaussian kernel for a finite point source
    kernel = np.exp(-(temp_x**2 + temp_y**2)/2)  
    
    return(kernel/np.sum(kernel))



def path2displacement(path_vec):
    [h0, v0] = path_vec[0]
    h_delta = path_vec[:,0] - h0
    v_delta = path_vec[:,1] - v0
    return(np.sqrt(h_delta ** 2 + v_delta ** 2))


def bin_rays(i, path_point, global_beta, new_beta_boundaries, s_img, delta_beta):
    rel_beta_pixels = np.nonzero(((path_point[0] - global_beta[:, 0])**2 + (path_point[1] - global_beta[:, 1])**2)
                                < new_beta_boundaries**2)[0]
    # Magnification assignment
    if len(rel_beta_pixels) > 0:
        lc_value = len(rel_beta_pixels)*s_img/(np.pi*delta_beta**2)
    else:  
        lc_value = 0            
    return(lc_value)
