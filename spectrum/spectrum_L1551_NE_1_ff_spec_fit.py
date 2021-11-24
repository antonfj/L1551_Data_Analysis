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
init_guess = [2.0, 10.0, 1.0e-3, 0.01]
popt, pcov = scipy.optimize.curve_fit(combined_power_law_fit, freq, log_flux, init_guess)
print(popt)
alpha_high, K_1, K_2, K_3 = popt
print("High frequency spec. ind.: ", alpha_high)

############## 6. Plot all the spectra in graphs ############################### 
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
#fig = plt.figure(5, figsize=[9.6,6.3])
fig = plt.figure(5)
ax = fig.add_subplot(111)
plt.errorbar(log_freq, log_flux, yerr=log_flux_err, fmt='bo')

# Plot spectral fit
# Make array of frequencies to plot spectral fit
freq_samples = np.arange(4.0, 400.0, 0.1)
spectral_fit = combined_power_law_fit(freq_samples, alpha_high, K_1, K_2, K_3)

# Plot combined spectral fit over all frequencies
plt.plot(np.log10(freq_samples), spectral_fit, 'k--')

# Plot low freq. spectral fit
plt.plot(np.log10(freq_samples), np.log10((K_1*freq_samples**2)*(1 - np.exp(K_2*freq_samples**-2.1))), 'r-')

# Plot high freq. spectral fit
plt.plot(np.log10(freq_samples), alpha_high*np.log10(freq_samples) + np.log10(K_3), 'b-')

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
plt.savefig('spectrum_L1551_NE_1_ff_spec_fit.png')

plt.show()
