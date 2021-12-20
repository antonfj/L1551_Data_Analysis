import numpy as np
import matplotlib
import os
import aplpy

# Images to use in the final image
c_band_image = 'L1551_NE_C_Band_VLA_initial_r_+0.5.pbcor.cropped.fits'

final_image = 'L1551_NE_C_Band_VLA.png'

fig = aplpy.FITSFigure(c_band_image)
fig.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig.show_colorscale(stretch='linear')
fig.show_colorscale(vmin=-36e-6, vmax=0.33e-3, cmap='jet')

# Add colourbar to image
fig.add_colorbar()
fig.colorbar.set_axis_label_text(r'Flux (Jy/beam)')

# Overplot Ku Band contours
sigma=12e-6
fig.show_contour(c_band_image, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	10*sigma, 20*sigma, 30*sigma, 40*sigma, 50*sigma],
	colors='white', overlap=True)

# Adds synthesis beam in bottom left corner
fig.add_beam()
fig.beam.set_color('white')
fig.beam.set_pad(1.0)

#Set frame colour to black
fig.frame.set_color('white')
#Set tick colour to black
fig.ticks.set_color('white')

fig.save(final_image, dpi=500)

