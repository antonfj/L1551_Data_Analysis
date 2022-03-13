import numpy as np
import matplotlib
import os
import aplpy

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
simulated_image = 'L1551_IRS_5_C_Band_5_GHz_simulated.fits'
observed_image = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.fits'

final_image = 'L1551_IRS_5_C_Band_simulated_image_comparison.eps'

fig = matplotlib.pyplot.figure(figsize=[8.,3.5])

titlefontsize=12
labelfontsize=8
cmap='jet'
framecolor='white'
tickcolor='white'
labelcolor='white'
beamcolor='white'

###################################################
############### Observed image #################### 
###################################################
fig1 = aplpy.FITSFigure(observed_image, figure=fig, subplot=[0.13, 0.1, 0.35, 0.8])
fig1.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)

fig1.show_colorscale(stretch='linear')
fig1.show_colorscale(vmin=-2e-5, vmax=3.2e-4, cmap=cmap)

# Add colourbar to image
#fig.add_colorbar()
#fig.colorbar.set_axis_label_text(r'Flux (Jy/beam)')

# Set title
fig1.set_title('Observation', size=titlefontsize)

# Adds synthesis beam in bottom left corner
#fig.add_beam(major=0.0055, minor=0.0055, angle=90.0)
fig1.add_beam()
fig1.beam.set_color(beamcolor)
fig1.beam.set_pad(1.0)

# Set font size of labels
fig1.axis_labels.set_font(size=labelfontsize)
fig1.tick_labels.set_font(size=labelfontsize)

#Set frame colour
fig1.frame.set_color('white')
#Set tick colour
fig1.ticks.set_color('white')

###################################################
############### Simulated image ################### 
###################################################
fig2 = aplpy.FITSFigure(simulated_image, figure=fig, subplot=[0.47, 0.1, 0.45, 0.8])
fig2.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)
fig2.show_colorscale(stretch='linear')
fig2.show_colorscale(vmin=-2e-5, vmax=3.2e-4, cmap=cmap)

# Add colourbar to image
fig2.add_colorbar()
fig2.colorbar.set_axis_label_text(r'Flux (Jy/beam)')

# Set title
fig2.set_title('Model', size=titlefontsize)

# Adds synthesis beam in bottom left corner
#fig.add_beam(major=0.0055, minor=0.0055, angle=90.0)
fig2.add_beam()
fig2.beam.set_color(beamcolor)
fig2.beam.set_pad(1.0)

# Set font size of labels
fig2.axis_labels.set_font(size=labelfontsize)
fig2.tick_labels.set_font(size=labelfontsize)

# Remove Dec labels
fig2.axis_labels.hide_y()
fig2.tick_labels.hide_y()

#Set frame colour
fig2.frame.set_color('white')
#Set tick colour
fig2.ticks.set_color('white')

# Overplot C Band contours
#sigma=11e-6
#fig.show_contour(c_band_image, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
#	6*sigma, 7*sigma, 8*sigma, 9*sigma, 10*sigma, 15*sigma, 20*sigma, 25*sigma],
#	colors='white', overlap=True, slices=[0])

# Coordinates of sources
#north_coord = np.array([67.8923542,18.1346203])
#south_coord = np.array([67.8923708,18.1345203])

# Jet axis of sources
#jet_axis_north = 67.0 * (np.pi/180.0)
#jet_axis_south = 55.0 * (np.pi/180.0)

# Plot line along jet axes of both sources
#ra_length = 2.0 / (60*60)	# Set length of line along ra axis
#dec_length_north = ra_length / np.tan(jet_axis_north)
#dec_length_south = ra_length / np.tan(jet_axis_south)

#jet_coords = np.array([[ [north_coord[0] - ra_length/2, north_coord[0] + ra_length/2],
#	[north_coord[1] - dec_length_north/2, north_coord[1] + dec_length_north/2] ],
#	[ [south_coord[0] - ra_length/2, south_coord[0] + ra_length/2],
#	[south_coord[1] - dec_length_south/2, south_coord[1] + dec_length_south/2] ] ])
#fig.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'])

# Plot marker at position of both sources
#fig.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)


fig.savefig(final_image, dpi=500)

matplotlib.pyplot.show()