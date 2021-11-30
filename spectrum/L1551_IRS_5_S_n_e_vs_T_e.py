import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize

from uncertainties import ufloat
from uncertainties.umath import *  # Imports sin(), etc.

T_e_values = np.arange(0.1e4, 5e4, 0.1e4)  # Range of electron temperatures to test
print(T_e_values)

theta_values = np.arange(0.01, 0.5, 0.01)
print(theta_values)

#theta_maj = ufloat(0.396, 0.023)	# Major axis of source in arcsec
#theta_min = ufloat(0.057, 0.022) 	# Minor axis of source in arcsec
theta_maj = 0.202               	# Major axis of source in arcsec
theta_min = 0.180                       # Minor axis of source in arcsec

D = 0.14		# Distance to source in kpc

Omega_s = (math.pi*theta_maj*theta_min)/(4*math.log(2)) 	# Solid angle size of emission
print("Omega_s: ", Omega_s)

# Function for finding reduced chi-squared values
def reduced_chisquared(observed, expected, errors, no_parameters):
	dof = len(observed) - no_parameters
	sum_square_deviations = 0.0

	for i in range(len(observed)):
		square_deviation = (observed[i] - expected[i])**2 / errors[i]**2
		sum_square_deviations += square_deviation
	
	red_chi_squared = sum_square_deviations/dof
	return red_chi_squared
	
# Combined power law spectrum
def combined_power_law_fit(freq, alpha_high, K_2, K_3,):
        log_flux = np.log10((7.21586e-4*T_e*Omega_s*freq**2)*(1 - np.exp(-K_2*T_e**-1.35*freq**-2.1)) + K_3*freq**alpha_high)
        return log_flux

# Flux values
#flux = np.array([0.11, 0.51, 0.90, 1.19, 1.27, 1.63, 12.15, 34.07, 85.10, 271., 497.])
#flux_err = np.array([0.00, 0.10, 0.03, 0.03, 0.05, 0.04, 0.10, 0.76, 1.70, 17., 35.])
#flux = np.array([0.11, 0.51, 0.90, 1.19, 1.27, 1.63, 10.57, 13.63, 34.07, 85.10, 271., 497.])
#flux_err = np.array([0.00, 0.10, 0.03, 0.03, 0.05, 0.04, 0.20, 0.45, 0.76, 1.70, 17., 35.])
flux = np.array([0.11, 0.51, 0.90, 1.19, 1.27, 1.63, 34.07, 85.10, 271., 497.])
flux_err = np.array([0.00, 0.10, 0.03, 0.03, 0.05, 0.04, 0.76, 1.70, 17., 35.])

# Flux values from Rodriguez et al. (1998)
flux_1998 = np.array([0.70, 0.78, 1.2, 2.0, 7.4, 45.0])
flux_err_1998 = np.array([0.20, 0.02, 0.10, 0.2, 0.5, 6.0])

# Flux values from Rodriguez et al. (2003)
flux_2003 = np.array([0.51])
flux_err_2003 = np.array([0.02])

# Add in quadrature 10% absolute flux calibration error
flux_err = np.sqrt(flux_err**2 + (0.1*flux)**2)
print("Errors: ", flux_err)
flux_err_1998 = np.sqrt(flux_err_1998**2 + (0.1*flux_1998)**2)
print("Errors (1998): ", flux_err_1998)
flux_err_2003 = np.sqrt(flux_err_2003**2 + (0.1*flux_2003)**2)
print("Errors (1998): ", flux_err_2003)

# Logs of flux densities and errors for the log-log graph
log_flux = np.log10(flux)
log_flux_err = (flux_err) / (flux * np.log(10))
log_flux_1998 = np.log10(flux_1998)
log_flux_err_1998 = (flux_err_1998) / (flux_1998 * np.log(10))
log_flux_2003 = np.log10(flux_2003)
log_flux_err_2003 = (flux_err_2003) / (flux_2003 * np.log(10))

# Frequencies that flux was measured at (X data)
#freq = np.array([5.0, 10.0, 13.5, 17.5, 20.0, 24.0, 43.0, 93.0, 153.0, 225.0, 336.0])
#freq = np.array([5.0, 10.0, 13.5, 17.5, 20.0, 24.0, 41.0, 45.0, 93.0, 153.0, 225.0, 336.0])
freq = np.array([5.0, 10.0, 13.5, 17.5, 20.0, 24.0, 93.0, 153.0, 225.0, 336.0])

# Frequencies of 1998 data
wavelengths_1998 = np.array([180.0, 36.0, 20.0, 13.0, 7.0, 2.7])
freq_1998 = 3e8/(wavelengths_1998*1e6)
print("Frequencies (1998): ", freq_1998)

# Frequencies of 2003 data
wavelengths_2003 = np.array([35.0])
freq_2003 = 3e8/(wavelengths_2003*1e6)
print("Frequencies (2003): ", freq_2003)

# Log of frequencies for the log-log graph
log_freq = np.log10(freq)
log_freq_1998 = np.log10(freq_1998)
log_freq_2003 = np.log10(freq_2003)

############## 2. Find electron density for each T_e value #################
################################################################################
n_e_values = np.zeros(len(T_e_values))
print(n_e_values)
for i in range(len(T_e_values)):
    T_e = T_e_values[i]
    init_guess = [2.0, 1e4, 0.01]
    popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq, log_flux, init_guess)
    perr = np.sqrt(np.diag(pcov))

    #K_2 = ufloat(popt[1], perr[1])
    K_2 = popt[1]

    n_e = sqrt( 1710.0 * K_2 * D**(-1) * theta_min**(-1) )
    n_e_values[i] = n_e
    #print("Electron density (cm^-3): ", n_e)

    freq_c = (K_2*T_e**-1.35)**(1/2.1)				# Turnover freq. in GHz
    #print("Turnover frequency (GHz): ", freq_c)
    

############## 3. Plot n_e vs T_e ##################
#####################################################################

#Makes ticks point inwards
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
# Makes all text non-italic
matplotlib.rcParams.update({'mathtext.default': 'regular'})
# Sets font size of all text on graph
matplotlib.rcParams.update({'font.size': 16})
plt.rc('legend', fontsize=16)

# Plot data values
fig1 = plt.figure(1, figsize=(15, 9))
ax = fig1.add_subplot(111)
plt.errorbar(T_e_values, n_e_values, fmt='bo')

plt.xlabel('T_e (K)')
plt.ylabel(r'$n_{e}\ (cm^{-3})$')

# Set axes
plt.axis([0.0, 50000.0, 0.0, 150000.0])

plt.tight_layout()				# Make everything fit in window
plt.savefig('L1551_IRS_5_S_ff_spec_fit_n_e_vs_T_e.png')

############## 4. Find electron density for each theta value #################
################################################################################
T_e = 2e4
n_e_values = np.zeros(len(theta_values))
print(n_e_values)
for i in range(len(theta_values)):
    theta_maj = theta_values[i]
    theta_min = theta_values[i]
 
    Omega_s = (math.pi*theta_maj*theta_min)/(4*math.log(2)) 	# Solid angle size of emission 
    
    init_guess = [2.0, 1e4, 0.01]
    popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq, log_flux, init_guess)
    perr = np.sqrt(np.diag(pcov))

    #K_2 = ufloat(popt[1], perr[1])
    K_2 = popt[1]

    n_e = sqrt( 1710.0 * K_2 * D**(-1) * theta_min**(-1) )
    n_e_values[i] = n_e
    #print("Electron density (cm^-3): ", n_e)

    freq_c = (K_2*T_e**-1.35)**(1/2.1)				# Turnover freq. in GHz
    #print("Turnover frequency (GHz): ", freq_c)
    
############## 5. Plot n_e vs theta ##################
#####################################################################

#Makes ticks point inwards
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
# Makes all text non-italic
matplotlib.rcParams.update({'mathtext.default': 'regular'})
# Sets font size of all text on graph
matplotlib.rcParams.update({'font.size': 16})
plt.rc('legend', fontsize=16)

# Plot data values
fig2 = plt.figure(2, figsize=(15, 9))
ax = fig2.add_subplot(111)
plt.errorbar(theta_values, n_e_values, fmt='bo')

plt.xlabel(r'$\theta$ (arcsec)')
plt.ylabel(r'$n_{e}\ (cm^{-3})$')

# Set axes
plt.axis([0.0, 0.5, 0.0, 200000.0])

plt.tight_layout()				# Make everything fit in window
plt.savefig('L1551_IRS_5_S_ff_spec_fit_n_e_vs_theta.png')
plt.show()
