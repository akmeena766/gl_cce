import numpy as np

import input
from constants import *


def light_curve(Dos, source_plane_FOV, source_size, tEin_scale, mu_source, shear_parallel):
    # pixel size in meters in source plane
    pixel2meter = Dos*source_plane_FOV*tEin_scale*ANGLE_ARCSEC/input.beta_res
    
    # half width of Gaussian luminosity distribution for point source, in pixels. We define the source size as 3 sigma.
    half_width = int(source_size/pixel2meter)

    # Initial and final position in pixels, x-rows, y-columns
    if shear_parallel:
        x0 = half_width + 1  
        xf = input.beta_res - (half_width + 2)
        y0 = int(input.beta_res/2)
        yf = y0
    else:
        x0 = int(input.beta_res/2)
        xf = x0
        y0 = half_width + 1  
        yf = input.beta_res - (half_width + 2)

    # Source velocity (in Einstein_angle/second)
    vt_sEin = input.vt/Dos/ANGLE_ARCSEC/tEin_scale		

    # Source pixel crossing time in years
    pixel2year = ((source_plane_FOV/input.beta_res)/vt_sEin)/year2sec

    # Source kernel
    WPK_kernel = gaussian_kernel(half_width)
    
    time_steps, light_curve = path_conv(mu_source, x0, y0, xf, yf, WPK_kernel, pixel2year)

    # to make sure light curves with different half-size match on the same axis
    time_steps += pixel2year*(3*half_width + 1)  

    return(time_steps, light_curve)
    


def gaussian_kernel(half_width):
    if half_width == 0:
        # point source
        kernel = np.zeros((3, 3))
        kernel[1, 1] = 1
    else:
        # size larger than pixel
        temp_x, temp_y = np.meshgrid(np.arange(-half_width, half_width + 1, 1),
                                     np.arange(-half_width, half_width + 1, 1))
        
        # a Gaussian kernel for a finite point source
        kernel = np.exp(-(temp_x**2 + temp_y**2)/2)  

    return(kernel/np.sum(kernel))
    
    
    
def path_conv(mu, x0, y0, xf, yf, kernel, pixel2year):
    # Number of grid points along x and y axis in magnification map
    nx, ny = mu.shape[0], mu.shape[1]

    # Number of grid points along x and y axis in kernel
    nkx, nky = kernel.shape[0], kernel.shape[1]

    dx=int(nkx // 2)
    dy=int(nky // 2)

    if dx+1+max(x0,xf)>nx or min(x0,xf)-dx<0 or dy+1+max(y0,yf)>ny or min(y0,yf)-dy<0:
        print('Kernel too large for chosen path')
        return 0
    
    # Normalization factor for each kernel
    norm=np.sum(kernel) 
    
    # 
    xpath, ypath = path_finder(x0, y0, xf, yf)

    convolved_path=[]
    for t in range(xpath.shape[0]):
        temp_arr = mu[ypath[t]-dy:ypath[t]+dy+1, xpath[t]-dx:xpath[t]+dx+1]
  
        temp_conv = np.sum(temp_arr*kernel)/norm
  
        convolved_path = np.append(convolved_path, temp_conv)

    time_path = pixel2year*np.arange(0, convolved_path.shape[0])
    return time_path,convolved_path
    


def path_finder(x0, y0, xf, yf):
    # If vertical path is chosen
    if x0 == xf:
        y = np.arange(y0, yf + 1, 1, dtype=int)
        x = x0 * np.ones(y.shape[0], dtype=int)
        return(x,y)

    # If parallel path is chosen
    if y0 == yf:
        x = np.arange(x0, xf + 1, 1, dtype=int)
        y = y0 * np.ones(x.shape[0], dtype=int)
        return(x,y)
    
    x=np.array([x0])
    y=np.array([y0])
    
    slope=(yf-y0)/(xf-x0)
    intercept=(xf*y0-x0*yf)/(xf-x0)
    
    while x[-1] != xf or y[-1] != yf:
        xtemp, ytemp = np.meshgrid(np.array([x[-1]-1, x[-1], x[-1]+1]), np.array([y[-1]-1, y[-1], y[-1]+1]))
        
        dtemp=np.sqrt((xtemp*slope - ytemp + intercept)**2/(slope**2 + 1)) + np.sqrt((xf - xtemp)**2 + (yf - ytemp)**2)
        idx = np.unravel_index(np.argmin(dtemp, axis=None), dtemp.shape)
        
        x = np.append(x, xtemp[idx])
        y = np.append(y, ytemp[idx])
    return(x,y)
