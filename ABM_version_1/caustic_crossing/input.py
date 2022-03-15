#### Cosmology parameters
OM_0 = 0.3					# Present Omega_matter
OB_0 = 0.0					# Present Omage_baryon
OL_0 = 0.7					# Present Omega_lambda

H0 = 70						# present Hubble constant


#### Gravitational lens systems redshifts
zol = 0.5					# Lens redshift
zos = 1.5					# Source redshift


#### Source parameters
source_size = 1.e12				# Source size (in meters)
vt = 1000*1.e3					# Source velocity (in meters/second)
light_curve_yr = 20 				# length of transverse motion of source (in years)
sub_light_curve_yr = 0.5			# Sub light curve size to calculate at a time
                               			# sub_light_curve_yr = 0 ===> One point at a time
                                		# sub_light_curve_yr > 0 ===> Bunch of points at once

#### Microlens population parameters
xbox = 0.20e-3					# Half box size in x-direction to distribute microlenses (in arcsec)
ybox = 10.00e-3					# Half box size in y-direction to distribute microlenses (in arcsec)
buff = 0.10e-3					# Buffer region (in arcsec)

stellar_density = 0.1				# Stellar density (in units of MASS_SUN/PC^2) 

lower_limit = 0.08				# Lower limit on the microlens mass
upper_limit = 1.50				# Upper limit on the microlens mass
alpha = 2.35					# Power law index (alpha = 2.35 for Salpeter IMF)
M_scale = 1.  					# Mass value (in solar mass) for scaling to make equations dimensionless


#### Microlens simulation parameters
user_seed = 12345				# Seed value for reproduction purpose

epsilon_ml_draw = 5.0	  			# The size increment factor of FOV in the image plane to draw microlenses
epsilon_ml_shoot = 3.0	  			# The size increment factor of FOV in the image plane to shoot the rays/pixels

method = 'ABM'  				# Method for calculation: 'IRS' or 'IPM' or 'ABM'

max_memory = 0.2 	 			# Maximum memory usage for arrays in intermediate calculations (in GB)


#### IRS settings
beta_res = 2000      				# Resolution in pixels per axis for plotting the source FOV (only for IRS and IPM)
rays_per_pixel = 10  			        # Number of rays per image plane pixels


#### ABM settings
n_pixels = 12 					# Number of image plane pixels per dimension in each division; could change to increase efficiency
eta_ratio = 0.7
n_steps = 1000    				# Number of initial light curve time steps (Half the number of final time steps)
beta_0 = 4  					# Factor of initial search boundary in the source plane (in units of FOV)

point_based = True				# If the 
kernel_flag = True  				# Whether to bin rays in the source plane and apply a Gaussian profile before creating light curve


#### Save data flag
save_data = True				# Whether to save relevant data

#### Plot flags
plot_critical = True				# Only for IRS and IPM 
plot_caustics = True				# Only for IRS and IPM
plot_light_curve_parallel = False		# Parallel light curve (along x-axis)
plot_light_curve_vertical = True		# Vertical light curve (along y-axis)
