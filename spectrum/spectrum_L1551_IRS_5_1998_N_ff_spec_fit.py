import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize

# Function for finding reduced chi-squared values
def reduced_chisquared(observed, expected, errors, no_parameters):
	dof = len(observed) - no_parameters
#	print dof
	sum_square_deviations = 0.0

	for i in range(len(observed)):
		square_deviation = (observed[i] - expected[i])**2 / errors[i]**2
		sum_square_deviations += square_deviation
#		print square_deviation
#		print sum_square_deviations
	
	red_chi_squared = sum_square_deviations/dof
#	print red_chi_squared
	return red_chi_squared

# Combined power law spectrum
def combined_power_law_fit(freq, alpha_high, K_1, K_2, K_3):
	log_flux = np.log10((K_1*freq**2)*(1 - np.exp(K_2*freq**-2.1)) + K_3*freq**alpha_high)
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

############## 2. Find spectral index in each part of spectrum #################
################################################################################
init_guess = [2.0, 10.0, 1.0e-3, 0.01]
popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq_1998, log_flux_1998, init_guess)
print(popt)
alpha_high, K_1, K_2, K_3 = popt
#print("Low frequency spec. ind.: ", alpha_low)
print("High frequency spec. ind.: ", alpha_high)


############## 6. Plot all the spectra in graphs ############################### 
################################################################################

# Make major tick labels in log-log graph corresponding to frequencies listed in following array
new_x_tick_labels=np.array([1, 5, 10, 50, 100, 500])
# Minor tick labels every 1 GHz from 1 - 10 GHz, every 10 GHz from 10 - 100 GHz, and every 100 GHz from 100 - 500 GHz
new_x_minor_tick_labels = np.concatenate((np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 500.0, 100.0)))
# Get log of frequencies for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_x_tick_locations=np.log10(new_x_tick_labels*(1.0))
new_x_minor_tick_locations=np.log10(new_x_minor_tick_labels*(1.0))

# Make major tick labels in log-log graph corresponding to flux densities listed in following array
new_y_tick_labels=np.array([0.01, 0.1, 1, 10, 50, 100, 500, 1000])
# Minor tick labels every 10 uJy from 10 - 100 uJy and every 100 uJy from 100 - 1500 uJy
new_y_minor_tick_labels=np.concatenate((np.arange(0.01, 0.1, 0.01), np.arange(0.1, 1.0, 0.1), np.arange(1.0, 10.0, 1.0), np.arange(10.0, 100.0, 10.0), np.arange(100.0, 1000.0, 100.0)))
# Get log of flux densities for tick labels so they can be applied to get corresponding tick locations in the log-log graph
new_y_tick_locations=np.log10(new_y_tick_labels)
new_y_minor_tick_locations=np.log10(new_y_minor_tick_labels)


############## 8. Plot the spectrum with all fits  ##################
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
fig = plt.figure(5, figsize=(15, 9))
ax = fig.add_subplot(111)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo', label='2020')
# Plot 2003 data
plt.errorbar(log_freq_2003, log_flux_2003, yerr=log_flux_err_2003, fmt='gs', label='2003')
# Plot 1998 data
plt.errorbar(log_freq_1998, log_flux_1998, yerr=log_flux_err_1998, fmt='r^', label='1998')

# Plot spectral fit
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(1.0, 400.0, 0.1)
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_1, K_2, K_3)

# Plot combined spectral fit over all frequencies
plt.plot(np.log10(freq_samples), spectral_fit, 'k--', label='Combined')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), np.log10((K_1*freq_samples**2)*(1 - np.exp(K_2*freq_samples**-2.1))), color='green', label='Free-Free')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), color='brown', label='Dust')

# Set axes for the log-log graph by getting log of upper and lower limits for axes
plt.axis(np.log10([1.0, 500.0, 0.01, 1000.0]))

# Apply major and minor ticks for x and y axes
ax.set_xticks(new_x_tick_locations)
ax.set_xticklabels(new_x_tick_labels)
ax.set_xticks(new_x_minor_tick_locations, minor=True)
ax.set_yticks(new_y_tick_locations)
ax.set_yticklabels(new_y_tick_labels)
ax.set_yticks(new_y_minor_tick_locations, minor=True)

plt.xlabel(r'$\mathit{\nu}$ (GHz)')
plt.ylabel(r'$\mathit{S_{\nu}}\ (mJy$)')

# Add legend to plot
plt.legend(loc='upper left')

plt.tight_layout()				# Make everything fit in window
plt.savefig('spectrum_L1551_IRS_5_N_1998_ff_spec_fit_with_2020_data.png')

plt.show()
