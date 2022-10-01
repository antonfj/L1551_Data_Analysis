"""
Basic docstring explaining example
"""
from __future__ import print_function
#********************
#sf3dmodels libraries
#********************
from sf3dmodels.outflow import OutflowModel   #Model functions
import sf3dmodels.utils.units as u            #Units
import sf3dmodels.rt as rt                    #Writing functions for radiative transfer
import sf3dmodels.Plot_model as Pm            #Plotting model
import sf3dmodels.Model as Model              #Grid
from sf3dmodels.grid import Overlap           #Overlap submodels
#********************
#radmc3dpy libraries
#********************
from radmc3dPy.image import *
#********************
#Extra libraries
#********************
import numpy as np
import time
import subprocess
from matplotlib import cm

################################################
# Set parameters of outflow to be modelled
################################################

t0 = time.time()
#--------------------

dx_grid = 0.25*u.au
#*********************
#OUTFLOW 1
#*********************
tag = '_outflow1'

#---------------------
#GEOMETRIC PARAMETERS
#---------------------
pos_c = np.array([0*u.au, 0*u.au, 0])
axis = np.array([-1,0.7,-1.22]) 
# Min and max distances to calculate jet parameters
z_min = 10*u.au
z_max = 120*u.au

z0 = z_min               # Injection radius of jet

op_angle = 0.489         # Opening angle of jet in radians

#---------------------
#PHYSICAL PROPERTIES
#---------------------
#w0 = 1*u.au
w0 = z0 * np.tan(op_angle/2)
print("w0 (au): ", w0/u.au)
eps = 0.66
w = [w0, eps]

T0 = 1e4
qT = 0.0
temp = [T0, qT]

ionfrac = [1.0,0.0]

abund = [1e-4, 0]
gtd = 100

v0 = 200 * 1e3 #km/s

mass_loss_rate=4.7e-9     # Mass loss rate of star in M_sun yr^-1

#dens0 = 2e12
dens0 = (1.2e49/(w0**2 * v0)) * mass_loss_rate # Calculate initial density from mass loss rate
print("n_0: {:.2e}".format(dens0))
qn = -2*eps
dens = [dens0, qn]

#---------------------
#COMPUTING MODEL
#---------------------
Outf1 = OutflowModel(pos_c, axis, z_min, z_max, dx_grid) #Initializing Class with grid parameters
Outf1.reynolds86(w, dens, ionfrac, temp, v0, abund, gtd, z0=z0) #Invoking the outflow model from Reynolds et al. 1986

#---------------------
#WRITING FOR RADMC3D
#---------------------
#Using the Radmc3d class
prop1 = {'dens_e': Outf1.density_ion,
         'dens_ion': Outf1.density_ion,
         'temp_gas': Outf1.temperature}
radmc1 = rt.Radmc3d(Outf1.GRID)
radmc1.submodel(prop1, output='datatab'+tag+'.dat')
print('Output columns', radmc1.columns)

#-------
#TIMING
#-------
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

###################################################
# Model outflow with sf3dmodels
###################################################

t0 = time.time()

#********
#GRIDDING
#********
sizex = 200 * u.au
sizey = sizez = 200 * u.au
Nx = Ny = Nz = 800
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz], rt_code='radmc3d')

#********
#MERGING
#********
files = ['datatab_outflow1.dat']
#Old style
#outflows = BGG.overlap(GRID, submodels = data2merge, rho_min = 1e6)

#New style
columns = ['id', 'x', 'y', 'z', 'dens_e', 'dens_ion', 'temp_gas']
outflows = Overlap(GRID)
finalprop = outflows.fromfiles(columns, submodels = files, rt_code = 'radmc3d')

radmc = rt.Radmc3dDefaults(GRID)
radmc.freefree(finalprop)
#********
#TIMING
#********
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

#********
#PLOTTING
#********
density = finalprop['dens_e'] / 1e6 #dens. in cm^-3
temperature = finalprop['temp_gas']

weight = 100 * np.mean(density)

#-----------------
#Plot for DENSITY
#-----------------
Pm.scatter3D(GRID, density, weight, NRand = 4000, axisunit = u.au, colorscale = 'log', cmap = 'cool',
             colorlabel = r'${\rm log}_{10}(n [cm^{-3}])$', output = 'global_grid_dens.png', vmin = 5, show=False)

#--------------------
#Plot for TEMPERATURE
#--------------------
Pm.scatter3D(GRID, density, weight, colordim = temperature, NRand = 4000, axisunit = u.au, colorscale = 'log',
             cmap = 'brg', colorlabel = r'${\rm log}_{10}(T$ $[K])$', output = 'global_grid_temp.png', vmin = 2, show=False)

#####################################################################
# Calculate radiative transfer and simulate image with RADMC-3D
#####################################################################

t0 = time.time()

num_threads = "24"      # Set num. of threads for RADMC-3D to use

# Image at C Band (5 GHz)
frequency = 5e9         # Frequency to simulate image
sizeau = "270"          # Size in au of area to image
npix = "100"             # Size of image in pixels

# Image at C Band (5 GHz)
radmc3d_command = ["radmc3d", "image", "lambda", "60000", "setthreads", num_threads, "sizeau", sizeau, "npix", npix, "stokes"]
subprocess.run(radmc3d_command, shell=False, check=True, universal_newlines=True)
#makeImage(wav=60000.,sizeau=sizeau,npix=npix)

#********
#TIMING
#********
print ('-------------------------------------------------\n')
print ('Ellapsed time: %.3fs' % (time.time() - t0))
print ('-------------------------------------------------\n-------------------------------------------------\n')

C_Band_image = readImage("image.out")
# Convolve beam of e-MERLIN C Band image
C_Band_image = C_Band_image.imConv(dpc=140.,psfType='gauss',fwhm=[0.220,0.113],pa=43.831)
#plotImage(C_Band_image,log=False,cmap=cm.hot,bunit='jy/beam',dpc=140,au=False,arcsec=True)

#Create fits file of image at coordinates of L1551 IRS 5
C_Band_image.writeFits(fname="L1551_IRS_5_C_Band_5_GHz_simulated_mlr_4.7e-9_z0_10au_zmax_120au_res_0.25au_grid_0.25au_coarse.fits", dpc=140,coord="04h31m34.17s+18d08m04.28s",
                       bandwidthmhz=500.0, casa=False,stokes="I")

"""
# Same image convolved to X Band beam
C_Band_image = readImage("image.out")
C_Band_image = C_Band_image.imConv(dpc=140.,psfType='gauss',fwhm=[0.202,0.180],pa=-43.027)
plotImage(C_Band_image,log=False,cmap=cm.hot,bunit='jy/beam',dpc=140,au=False,arcsec=True)
"""