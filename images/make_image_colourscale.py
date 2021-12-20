import numpy as np
import matplotlib
import os
import aplpy

# Images to use in the final image
alma_image = 'L1551_IRS_5_ALMA_Band_4_shifted_to_VLA_epoch.fits'
vla_ku_band_image = 'L1551_IRS5_Ku_Band_VLA_r_0_5_sc1.cropped.fits'

final_image = 'L1551_IRS_5_ALMA+Ku_Band.png'

fig = aplpy.FITSFigure(alma_image)
fig.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig.show_colorscale(stretch='linear')
fig.show_colorscale(vmin=-2e-4, vmax=1.3e-2, cmap='jet')

# Add colourbar to image
fig.add_colorbar()
fig.colorbar.set_axis_label_text(r'Flux (Jy/beam)')

"""
# Add labels for bow shock, counter-jet and star
fig.add_label(0.25, 0.6, 'Counter-Jet', relative=True, color='black')
fig.add_label(0.65, 0.1, 'Bow Shock', relative=True, color='black')
fig.add_label(0.58, 0.6, 'DG Tau A', relative=True, color='black')
"""

# Overplot Ku Band contours
sigma=4.3e-6
fig.show_contour(vla_ku_band_image, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	6*sigma, 7*sigma, 8*sigma, 9*sigma, 10*sigma, 20*sigma, 40*sigma, 60*sigma,
	80*sigma, 100*sigma, 120*sigma, 140*sigma, 160*sigma, 180*sigma], colors='white', overlap=True)

# Coordinates of sources
north_coord = np.array([67.8923542,18.1346203])
south_coord = np.array([67.8923708,18.1345203])

# Jet axis of sources
jet_axis_north = 67.0 * (np.pi/180.0)
jet_axis_south = 55.0 * (np.pi/180.0)

# Plot line along jet axis
ra_length = 2.0 / (60*60)	# Set length of line along ra axis
dec_length_north = ra_length / np.tan(jet_axis_north)
dec_length_south = ra_length / np.tan(jet_axis_south)

jet_coords = np.array([[ [north_coord[0] - ra_length/2, north_coord[0] + ra_length/2],
	[north_coord[1] - dec_length_north/2, north_coord[1] + dec_length_north/2] ],
	[ [south_coord[0] - ra_length/2, south_coord[0] + ra_length/2],
	[south_coord[1] - dec_length_south/2, south_coord[1] + dec_length_south/2] ] ])
fig.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'])

# Plot marker at position of northern source
fig.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig.add_beam()
fig.beam.set_color('white')
fig.beam.set_pad(1.0)

#Set frame colour to black
fig.frame.set_color('white')
#Set tick colour to black
fig.ticks.set_color('white')

fig.save(final_image, dpi=500)

