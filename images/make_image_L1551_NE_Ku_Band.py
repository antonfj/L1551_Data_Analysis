import numpy as np
import matplotlib
import os
import aplpy

# Images to use in the final image
alma_image = 'L1551_NE_ALMA_Band_7.cont.I.pbcor.fits'
vla_ku_band_image = 'L1551_NE_Ku_Band_VLA_initial_r_+0.5.fits'

final_image = 'L1551_NE_ALMA_Band_7+VLA_Ku_Band.png'

fig = aplpy.FITSFigure(alma_image)
fig.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig.show_colorscale(stretch='linear')
fig.show_colorscale(vmin=-24e-4, vmax=0.2, cmap='jet')

# Add colourbar to image
fig.add_colorbar()
fig.colorbar.set_axis_label_text(r'Flux (Jy/beam)')

# Overplot Ku Band contours
sigma=4.3e-6
fig.show_contour(vla_ku_band_image, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	10*sigma, 20*sigma, 40*sigma, 60*sigma, 80*sigma], colors='white', overlap=True)

# Coordinates of sources
# Took peak of Ku Band emission as position of source
A_coord = np.array([67.93547224,18.14202352])
B_coord = np.array([67.93533817,18.14209124])

# Jet axis of sources
jet_axis_A = 228.0 * (np.pi/180.0)
jet_axis_B = 238.0 * (np.pi/180.0)

# Plot line along jet axis
ra_length = 2.0 / (60*60)	# Set length of line along ra axis
dec_length_A = ra_length / np.tan(jet_axis_A)
dec_length_B = ra_length / np.tan(jet_axis_B)

jet_coords = np.array([[ [A_coord[0] - ra_length/2, A_coord[0] + ra_length/2],
	[A_coord[1] - dec_length_A/2, A_coord[1] + dec_length_A/2] ],
	[ [B_coord[0] - ra_length/2, B_coord[0] + ra_length/2],
	[B_coord[1] - dec_length_B/2, B_coord[1] + dec_length_B/2] ] ])
fig.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'])

# Plot marker at positions of both sources
fig.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig.add_beam()
fig.beam.set_color('white')
fig.beam.set_pad(1.0)

#Set frame colour to black
fig.frame.set_color('white')
#Set tick colour to black
fig.ticks.set_color('white')

fig.save(final_image, dpi=500)

