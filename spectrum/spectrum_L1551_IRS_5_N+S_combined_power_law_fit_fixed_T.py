import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize

from uncertainties import ufloat
from uncertainties.umath import *  # Imports sin(), etc.

# Parameters of jet
T_e = 1e4                               # Fixed value of electron temperature
r_0 = 8                            # Radius of base of jet in au
w_0 = 1.4                             # Width at base of jet in au
i = 45. * (2*math.pi) / (360)  # Inclination angle of jet in radians
print("Inclination angle (radians): ", i)

#theta_maj = ufloat(0.396, 0.023)	# Major axis of source in arcsec
#theta_min = ufloat(0.057, 0.022) 	# Minor axis of source in arcsec
theta_maj = 0.454               	# Major axis of source in arcsec
theta_min = 0.195                       # Minor axis of source in arcsec

D = 140.                            # Distance to source in pc

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
def combined_power_law_fit(freq, alpha_high, K_1, K_3):
        log_flux = np.log10(K_1*freq**0.25 + K_3*freq**alpha_high)
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
new_y_tick_labels=np.array([0.01, 0.1, 1, 10, 100, 1000])
# Minor tick labels every 10 uJy from 10 - 100 uJy and every 100 uJy from 100 - 1500 uJy
new_y_minor_tick_labels=np.concatenate((np.arange(0.01,0.1,0.01),np.arange(0.1,1.0,0.1), np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 1000.0, 100.0)))
# Get log of flux densities for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_y_tick_locations=np.log10(new_y_tick_labels)
new_y_minor_tick_locations=np.log10(new_y_minor_tick_labels)

# Make y labels strings to not have decimal points
new_y_tick_labels=np.array(['0.01', '0.1',  '1', '10', '100', '1000'])

######################################################################################
############################### N source #############################################
######################################################################################

# Flux values
flux = np.array([0.11, 0.47, 0.88, 1.15, 1.23, 1.60, 32.2, 78.6, 240., 418.])
flux_err = np.array([0.02, 0.03, 0.05, 0.06, 0.07, 0.09, 1.66, 4., 12.8, 22.])

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

# Frequencies that flux was measured at
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

############## 8. Plot the spectrum with fit ########################
#####################################################################

# Plot data values
fig = plt.figure(1, figsize=(7.2, 3.2))
ax1 = fig.add_subplot(121)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo', label='2020', markersize=markersize)
# Plot 2003 data
plt.errorbar(log_freq_2003, log_flux_2003, yerr=log_flux_err_2003, fmt='gs', label='2003', markersize=markersize)
# Plot 1998 data
plt.errorbar(log_freq_1998, log_flux_1998, yerr=log_flux_err_1998, fmt='r^', label='1998', markersize=markersize)

# Plot spectral fit
alpha_high, K_1, K_3 = popt
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(1.0, 500.0, 0.1)

# Plot combined spectral fit over all frequencies
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_1, K_3)
plt.plot(np.log10(freq_samples), spectral_fit, 'k--', label='Combined')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), 0.25*np.log10(freq_samples) + np.log10(K_1), color='green', label='Free-Free')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), color='darkorange', label='Dust')

# Apply major and minor ticks for x and y axes
ax1.set_xticks(new_x_tick_locations)
ax1.set_xticklabels(new_x_tick_labels)
ax1.set_xticks(new_x_minor_tick_locations, minor=True)
ax1.set_yticks(new_y_tick_locations)
ax1.set_yticklabels(new_y_tick_labels)
ax1.set_yticks(new_y_minor_tick_locations, minor=True)

plt.xlabel(r'$\mathit{\nu}$ (GHz)')
plt.ylabel(r'$\mathit{S_{\nu}}\ (mJy$)')
plt.title('L1551 IRS 5 N')

# Add legend to plot
plt.legend(loc='upper left', borderaxespad=0.7)

plt.tight_layout()				# Make everything fit in window


######################################################################################
############################### S source #############################################
######################################################################################
# Flux values
flux = np.array([1.16, 1.75, 1.97, 2.05, 2.12, 2.10, 18.4, 42.8, 125., 210.])
flux_err = np.array([0.07, 0.10, 0.10, 0.11, 0.12, 0.12, 1.0, 2.2, 7., 12.])

# Flux values from Rodriguez et al. (1998)
flux_1998 = np.array([0.60, 0.70, 0.9, 1.5, 4.8, 23.0])
flux_err_1998 = np.array([0.20, 0.02, 0.10, 0.2, 0.5, 6.0])

# Flux values from Rodriguez et al. (2003)
flux_2003 = np.array([1.14])
flux_err_2003 = np.array([0.03])

# Add in quadrature 10% absolute flux calibration error
#flux_err = np.sqrt(flux_err**2 + (0.05*flux)**2)
#print("Errors: ", flux_err)

# Logs of flux densities and errors for the log-log graph
log_flux = np.log10(flux)
log_flux_err = (flux_err) / (flux * np.log(10))
log_flux_1998 = np.log10(flux_1998)
log_flux_err_1998 = (flux_err_1998) / (flux_1998 * np.log(10))
log_flux_2003 = np.log10(flux_2003)
log_flux_err_2003 = (flux_err_2003) / (flux_2003 * np.log(10))

# Frequencies that flux was measured at
freq = np.array([5.1, 10.0, 13.5, 17.5, 20.0, 24.0, 93.0, 153.0, 225.0, 336.0])

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

############## Plot the spectrum with fit ###########################
#####################################################################

# Plot data values
# Share axes with previous subplot
ax2 = fig.add_subplot(122, sharex=ax1, sharey=ax1)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo', label='2020', markersize=markersize)
# Plot 2003 data
plt.errorbar(log_freq_2003, log_flux_2003, yerr=log_flux_err_2003, fmt='gs', label='2003', markersize=markersize)
# Plot 1998 data
plt.errorbar(log_freq_1998, log_flux_1998, yerr=log_flux_err_1998, fmt='r^', label='1998', markersize=markersize)

# Plot spectral fit
alpha_high, K_1, K_3 = popt
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(1.0, 500.0, 0.1)

# Plot combined spectral fit over all frequencies
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_1, K_3)
plt.plot(np.log10(freq_samples), spectral_fit, 'k--', label='Combined')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), 0.25*np.log10(freq_samples) + np.log10(K_1), color='green', label='Free-Free')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), color='darkorange', label='Dust')

# Set axes for the log-log graph by getting log of upper and lower limits for axes
plt.axis(np.log10([1.0, 500.0, 0.01, 1000.0]))

# Make y axis labels invisible
ax2.tick_params(axis='y', labelleft=False)

plt.xlabel(r'$\mathit{\nu}$ (GHz)')
plt.title('L1551 IRS 5 S')

# Add legend to plot
plt.legend(loc='upper left', borderaxespad=0.7)     # Increase paddding between legend and axes

plt.tight_layout()				# Make everything fit in window

plt.savefig('spectrum_L1551_IRS_5_N+S_combined_power_law_fit_fixed_T_with_1998_data.pdf')

plt.show()
