import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize

from uncertainties import ufloat
from uncertainties.umath import *  # Imports sin(), etc.

T_e = 1e4                               # Fixed value of electron temperature
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
def combined_power_law_fit(freq, alpha_high, K_2, K_3):
        log_flux = np.log10((7.21586e-4*T_e*Omega_s*freq**2)*(1 - np.exp(-K_2*T_e**-1.35*freq**-2.1)) + K_3*freq**alpha_high)
        return log_flux

# Flux values
flux = np.array([0.32, 0.30, 0.55, 0.79, 1.04, 1.31, 372.0])
flux_err = np.array([0.06, 0.03, 0.02, 0.02, 0.04, 0.05, 15.0])

# Flux values from Rodriguez et al. (1998)
flux_1998 = np.array([0.60, 0.70, 0.9, 1.5, 4.8, 23.0])
flux_err_1998 = np.array([0.20, 0.02, 0.10, 0.2, 0.5, 6.0])


# Add in quadrature 10% absolute flux calibration error
flux_err = np.sqrt(flux_err**2 + (0.1*flux)**2)
print("Errors: ", flux_err)
flux_err_1998 = np.sqrt(flux_err_1998**2 + (0.1*flux_1998)**2)
print("Errors (1998): ", flux_err_1998)

# Logs of flux densities for the log-log graph
log_flux = np.log10(flux)
log_flux_err = (flux_err) / (flux * np.log(10))

# Frequencies that flux was measured at (X data)
freq = np.array([6.0, 10.0, 13.5, 17.5, 20.0, 24.0, 336.0])

# Log of frequencies for the log-log graph
log_freq = np.log10(freq)

############## 2. Find spectral index in each part of spectrum #################
################################################################################
init_guess = [2.0, 1e4, 0.01]
popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq, log_flux, init_guess)
perr = np.sqrt(np.diag(pcov))
perr = np.sqrt(np.diag(pcov))


alpha_high = ufloat(popt[0], perr[0])
K_2 = ufloat(popt[1], perr[1])
K_3 = ufloat(popt[2], perr[2])

print("K_2: ", K_2)
#print("Low frequency spec. ind.: ", alpha_low)
print("High frequency spec. ind.: ", alpha_high)

############## 3. Calculate electron density and ionized gas mass ###################
################################################################################
n_e = sqrt( 1710.0 * K_2 * D**(-1) * theta_min**(-1) )
print("Electron density (cm^-3): ", n_e)

freq_c = (K_2*T_e**-1.35)**(1/2.1)				# Turnover freq. in GHz
print("Turnover frequency (GHz): ", freq_c)

# Assume approximate spherical geometry for gas mass
r = 1.496e16 * D * (theta_maj/2)       # Calculate radius of gas
print("r (cm): ", r)

m_H = 1.673e-24                         # Mass of hydrogen atom in g
M_sun = 1.989e+33                       # Mass of sun in g
M_ion = (4/3) * math.pi * r**3 * m_H * n_e
print("M_ion (g): ", M_ion)
print("M_ion (solar masses): ", M_ion/M_sun)

############## 4. Plot all the spectra in graphs ############################### 
#############################################################################

# Make major tick labels in log-log graph corresponding to frequencies listed in following array
new_x_tick_labels=np.array([1.0, 5.0, 10.0, 50.0, 100.0, 500.0])
# Minor tick labels every 1 GHz from 1 - 10 GHz, every 10 GHz from 10 - 100 GHz, and every 100 GHz from 100 - 500 GHz
new_x_minor_tick_labels = np.concatenate((np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 500.0, 100.0)))
# Get log of frequencies for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_x_tick_locations=np.log10(new_x_tick_labels*(1.0))
new_x_minor_tick_locations=np.log10(new_x_minor_tick_labels*(1.0))

# Make major tick labels in log-log graph corresponding to flux densities listed in following array
new_y_tick_labels=np.array([0.1, 0.5, 1.0, 10.0, 50.0, 100.0, 500.0])
# Minor tick labels every 10 uJy from 10 - 100 uJy and every 100 uJy from 100 - 1500 uJy
new_y_minor_tick_labels=np.concatenate((np.arange(0.1,1.0,0.1), np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 500.0, 100.0)))
# Get log of flux densities for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_y_tick_locations=np.log10(new_y_tick_labels)
new_y_minor_tick_locations=np.log10(new_y_minor_tick_labels)


############## 5. Plot the spectrum with all fits  ##################
#####################################################################

#Makes ticks point inwards
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
# Makes all text non-italic
matplotlib.rcParams.update({'mathtext.default': 'regular'})
# Sets font size of all text on graph
matplotlib.rcParams.update({'font.size': 12})
plt.rc('legend', fontsize=12)

# Plot data values
fig = plt.figure(5, figsize=(7, 4))
ax = fig.add_subplot(111)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo')

# Plot spectral fit
alpha_high, K_2, K_3 = popt
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(1.0, 400.0, 0.1)

# Plot combined spectral fit over all frequencies
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_2, K_3)
plt.plot(np.log10(freq_samples), spectral_fit, 'k--', label='Combined')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), np.log10((7.21e-4*T_e*Omega_s*freq_samples**2)*(1 - np.exp(-K_2*T_e**-1.35*freq_samples**-2.1))), color='green', label='Free-Free')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), color='darkorange', label='Dust')

# Set axes for the log-log graph by getting log of upper and lower limits for axes
plt.axis(np.log10([1.0, 500.0, 0.1, 500.0]))

# Apply major and minor ticks for x and y axes
ax.set_xticks(new_x_tick_locations)
ax.set_xticklabels(new_x_tick_labels)
ax.set_xticks(new_x_minor_tick_locations, minor=True)
ax.set_yticks(new_y_tick_locations)
ax.set_yticklabels(new_y_tick_labels)
ax.set_yticks(new_y_minor_tick_locations, minor=True)

plt.xlabel(r'$\mathit{\nu}$ (GHz)')
plt.ylabel(r'$\mathit{S_{\nu}}\ (mJy$)')

plt.tight_layout()				# Make everything fit in window
plt.savefig('spectrum_L1551_NE_A_ff_spec_fit_fixed_T.png')

plt.show()
