import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize

from uncertainties import ufloat
from uncertainties.umath import *  # Imports sin(), etc.

# Parameters of jet
T_e = 1e4                           # Fixed value of electron temperature
r_0 = 10                            # Radius of base of jet in au

#w_0 = 1.4                           # Width at base of jet in au
i = 45. * (2*math.pi) / (360)  # Inclination angle of jet in radians
#print("Inclination angle (radians): ", i)

#theta_maj = ufloat(0.396, 0.023)	# Major axis of source in arcsec
#theta_min = ufloat(0.057, 0.022) 	# Minor axis of source in arcsec
theta_maj = 0.454               	# Major axis of source in arcsec
theta_min = 0.195                       # Minor axis of source in arcsec

D = 140.                            # Distance to source in pc

Omega_s = (math.pi*theta_maj*theta_min)/(4*math.log(2)) 	# Solid angle size of emission
#print("Omega_s: ", Omega_s)

alpha_low = 0.6                     # Spec. ind. of free-free emission

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
def combined_power_law_fit(freq, alpha_high, K_1, K_3):
        log_flux = np.log10(K_1*freq**alpha_low + K_3*freq**alpha_high)
        return log_flux

#####################################################################################
########################## Set figure parameters ####################################
#####################################################################################

fig = matplotlib.pyplot.figure(figsize=[7.2,3.2])

#Makes ticks point inwards
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
# Makes all text non-italic
matplotlib.rcParams.update({'mathtext.default': 'regular'})
# Sets font size of all text on graph
matplotlib.rcParams.update({'font.size': 8})
plt.rc('legend', fontsize=8)

markersize= 5        # set marker size

# Make major tick labels in log-log graph corresponding to frequencies listed in following array
new_x_tick_labels=np.array([1, 5, 10, 50, 100, 500])
# Minor tick labels every 1 GHz from 1 - 10 GHz, every 10 GHz from 10 - 100 GHz, and every 100 GHz from 100 - 500 GHz
new_x_minor_tick_labels = np.concatenate((np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 500.0, 100.0)))
# Get log of frequencies for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_x_tick_locations=np.log10(new_x_tick_labels*(1.0))
new_x_minor_tick_locations=np.log10(new_x_minor_tick_labels*(1.0))

# Make major tick labels in log-log graph corresponding to flux densities listed in following array
new_y_tick_labels=np.array([0.01, 0.1, 0.5, 1, 10, 50, 100, 500])
# Minor tick labels every 10 uJy from 10 - 100 uJy and every 100 uJy from 100 - 1500 uJy
new_y_minor_tick_labels=np.concatenate((np.arange(0.01,0.1,0.01),np.arange(0.1,1.0,0.1), np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 500.0, 100.0)))
# Get log of flux densities for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_y_tick_locations=np.log10(new_y_tick_labels)
new_y_minor_tick_locations=np.log10(new_y_minor_tick_labels)

# Make y labels strings to not have decimal points
new_y_tick_labels=np.array(['0.01', '0.1', '0.5', '1', '10', '50', '100', '500'])

######################################################################################
############################### A source #############################################
######################################################################################
print("L1551 NE A:")

# Flux values
flux = np.array([0.18, 0.31, 0.58, 0.81, 1.11, 1.36, 330.0])
flux_err = np.array([0.03, 0.04, 0.04, 0.05, 0.07, 0.08, 20.0])

# Flux values from Reipurth et al. (2002)
flux_2000 = np.array([0.39])
flux_err_2000 = np.array([0.01])


# Add in quadrature 10% absolute flux calibration error
flux_err = np.sqrt(flux_err**2 + (0.1*flux)**2)
#print("Errors: ", flux_err)
flux_err_2000 = np.sqrt(flux_err_2000**2 + (0.1*flux_2000)**2)
#print("Errors (2000): ", flux_err_2000)

# Logs of flux densities for the log-log graph
log_flux = np.log10(flux)
log_flux_err = (flux_err) / (flux * np.log(10))
log_flux_2000 = np.log10(flux_2000)
log_flux_err_2000 = (flux_err_2000) / (flux_2000 * np.log(10))

# Frequencies that flux was measured at (X data)
freq = np.array([6.0, 10.0, 13.5, 17.5, 20.0, 24.0, 336.0])

# Frequencies of 2000 data
wavelengths_2000 = np.array([36.0])
freq_2000 = 3e8/(wavelengths_2000*1e6)
#print("Frequencies (2000): ", freq_2000)

# Log of frequencies for the log-log graph
log_freq = np.log10(freq)
log_freq_2000 = np.log10(freq_2000)

############## 2. Find spectral index in each part of spectrum #################
################################################################################
init_guess = [2.0, 1, 0.01]
popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq, log_flux, init_guess)
perr = np.sqrt(np.diag(pcov))
perr = np.sqrt(np.diag(pcov))

#alpha_low = popt[0]
alpha_high = popt[0]
K_1 = popt[1]
K_3 = popt[2]
"""
# Give reduced chi-squared value (WARNING!!!: Reduced Chi-Square is not very accurate for non-linear fits)
log_expected_flux = combined_power_law_fit(freq, alpha_high, K_2, K_3)
print(log_expected_flux)
red_chi_squared = reduced_chisquared(log_flux, log_expected_flux, log_flux_err, 3)

print("Red. Chi-Squared: ", red_chi_squared)
"""
# Redefine fit parameters with errors
#alpha_low = ufloat(popt[0],perr[0])
alpha_high = ufloat(popt[0], perr[0])
K_1 = ufloat(popt[1], perr[1])
K_3 = ufloat(popt[2], perr[2])

print("K_1: ", K_1)
#print("Low frequency spec. ind.: ", alpha_low)
print("High frequency spec. ind.: ", alpha_high)

# Ionized Mass Loss Rate
# Using formula from Reynolds (1986)
jet_vel = 200 * (1e5)                   # Jet velocity in cm/s
op_angle = 30 * (2*math.pi) / (360)        # Opening angle of jet in radians

ionized_Mdot = 1.23e-18 * K_1**0.75 * op_angle**0.75 * jet_vel * D**1.5 * T_e**(-0.075) * (np.sin(i))**(-0.25)
print("Ionized Mass Loss Rate (M_sun/yr): ", ionized_Mdot)

############## Plot the spectrum with fit ###########################
#####################################################################

# Plot data values
fig = plt.figure(1, figsize=(7.2, 3.2))
ax1 = fig.add_subplot(121)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo', label='2020', markersize=markersize)

# Plot 2000 data
plt.errorbar(log_freq_2000, log_flux_2000, yerr=log_flux_err_2000, fmt='r^', label='2000', markersize=markersize)

# Plot spectral fit
alpha_high, K_1, K_3 = popt
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(1.0, 1000.0, 0.1)

# Plot combined spectral fit over all frequencies
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_1, K_3)
plt.plot(np.log10(freq_samples), spectral_fit, 'k--', label='Combined')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_low*np.log10(freq_samples) + np.log10(K_1), color='green', label='Free-Free')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), color='darkorange', label='Dust')

# Set axes for the log-log graph by getting log of upper and lower limits for axes
plt.axis(np.log10([1.0, 500.0, 0.01, 500.0]))

# Apply major and minor ticks for x and y axes
ax1.set_xticks(new_x_tick_locations)
ax1.set_xticklabels(new_x_tick_labels)
ax1.set_xticks(new_x_minor_tick_locations, minor=True)
ax1.set_yticks(new_y_tick_locations)
ax1.set_yticklabels(new_y_tick_labels)
ax1.set_yticks(new_y_minor_tick_locations, minor=True)

plt.xlabel(r'$\mathit{\nu}$ (GHz)')
plt.ylabel(r'$\mathit{S_{\nu}}\ (mJy$)')
plt.title('L1551 NE A')

# Add legend to plot
plt.legend(loc='upper left', borderaxespad=0.7)

plt.tight_layout()				# Make everything fit in window

######################################################################################
############################### B source #############################################
######################################################################################
print('\n')
print("L1551 NE B:")

# Flux values
flux = np.array([0.16, 0.15, 0.27, 0.33, 0.43, 0.51, 130.0])
flux_err = np.array([0.02, 0.03, 0.02, 0.02, 0.03, 0.04, 8.0])

# Flux values from Reipurth et al. (2002)
flux_2000 = np.array([0.27])
flux_err_2000 = np.array([0.01])

# Add in quadrature 10% absolute flux calibration error
flux_err = np.sqrt(flux_err**2 + (0.1*flux)**2)
flux_err_2000 = np.sqrt(flux_err_2000**2 + (0.1*flux_2000)**2)

# Logs of flux densities for the log-log graph
log_flux = np.log10(flux)
log_flux_err = (flux_err) / (flux * np.log(10))
log_flux_2000 = np.log10(flux_2000)
log_flux_err_2000 = (flux_err_2000) / (flux_2000 * np.log(10))

# Frequencies that flux was measured at (X data)
freq = np.array([6.0, 10.0, 13.5, 17.5, 20.0, 24.0, 336.0])

# Frequencies of 2000 data
wavelengths_2000 = np.array([36.0])
freq_2000 = 3e8/(wavelengths_2000*1e6)
#print("Frequencies (2000): ", freq_2000)

# Log of frequencies for the log-log graph
log_freq = np.log10(freq)
log_freq_2000 = np.log10(freq_2000)

############## 2. Find spectral index in each part of spectrum #################
################################################################################
init_guess = [2.0, 1, 0.01]
popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq, log_flux, init_guess)
perr = np.sqrt(np.diag(pcov))
perr = np.sqrt(np.diag(pcov))

#alpha_low = popt[0]
alpha_high = popt[0]
K_1 = popt[1]
K_3 = popt[2]
"""
# Give reduced chi-squared value (WARNING!!!: Reduced Chi-Square is not very accurate for non-linear fits)
log_expected_flux = combined_power_law_fit(freq, alpha_high, K_2, K_3)
print(log_expected_flux)
red_chi_squared = reduced_chisquared(log_flux, log_expected_flux, log_flux_err, 3)

print("Red. Chi-Squared: ", red_chi_squared)
"""
# Redefine fit parameters with errors
#alpha_low = ufloat(popt[0],perr[0])
alpha_high = ufloat(popt[0], perr[0])
K_1 = ufloat(popt[1], perr[1])
K_3 = ufloat(popt[2], perr[2])

print("K_1: ", K_1)
#print("Low frequency spec. ind.: ", alpha_low)
print("High frequency spec. ind.: ", alpha_high)

# Ionized Mass Loss Rate
# Using formula from Reynolds (1986)
jet_vel = 200 * (1e5)                   # Jet velocity in cm/s
op_angle = 30 * (2*math.pi) / (360)        # Opening angle of jet in radians

ionized_Mdot = 1.23e-18 * K_1**0.75 * op_angle**0.75 * jet_vel * D**1.5 * T_e**(-0.075) * (np.sin(i))**(-0.25)
print("Ionized Mass Loss Rate (M_sun/yr): ", ionized_Mdot)

############## Plot the spectrum with fit ###########################
#####################################################################

# Plot data values
# Share axes with previous subplot
ax2 = fig.add_subplot(122, sharex=ax1, sharey=ax1)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo', label='2020', markersize=markersize)

# Plot 2000 data
plt.errorbar(log_freq_2000, log_flux_2000, yerr=log_flux_err_2000, fmt='r^', label='2000', markersize=markersize)

# Plot spectral fit
alpha_high, K_1, K_3 = popt
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(1.0, 1000.0, 0.1)

# Plot combined spectral fit over all frequencies
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_1, K_3)
plt.plot(np.log10(freq_samples), spectral_fit, 'k--', label='Combined')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_low*np.log10(freq_samples) + np.log10(K_1), color='green', label='Free-Free')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), color='darkorange', label='Dust')

# Make y axis labels invisible
ax2.tick_params(axis='y', labelleft=False)

plt.xlabel(r'$\mathit{\nu}$ (GHz)')
plt.title('L1551 NE B')

# Add legend to plot
plt.legend(loc='upper left', borderaxespad=0.7)

plt.tight_layout()				# Make everything fit in window

plt.savefig('spectrum_L1551_NE_A+B_combined_power_law_fit_fixed_T_with_1998_data.pdf')

plt.show()
