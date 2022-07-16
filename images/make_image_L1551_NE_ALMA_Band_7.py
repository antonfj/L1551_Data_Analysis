import numpy as np
import matplotlib
import os
import aplpy
from astropy.io import fits

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
image = 'L1551_NE_ALMA_Band_7.cont.I.pbcor.fits'
image_mJy = 'L1551_NE_ALMA_Band_7.cont.I.pbcor.mJy.fits'

final_image = 'L1551_NE_ALMA_Band_7.pdf'

# Convert images from Jy to mJy
def convert_to_mJy(image, image_mJy):
    if not os.path.isfile(image_mJy):
        # Copy to new file
        os.system('cp ' + image + ' ' + image_mJy)

        # Opens FITS image
        hdulist = fits.open(image_mJy, mode='update')
        fitsdata = hdulist[0].data[:,0]

        # Converts to uJy
        fitsdata = fitsdata*1.e3

        hdulist[0].data[:,0] = fitsdata
        hdulist.flush()

convert_to_mJy(image,image_mJy)

titlefontsize=12
labelfontsize=8
cmap='jet'
framecolor='white'
tickcolor='white'
labelcolor='red'
beamcolor='white'

fig = aplpy.FITSFigure(image_mJy,figsize=[5.,3])
fig.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig.show_colorscale(stretch='linear')
fig.show_colorscale(vmin=-3, vmax=200, cmap='jet')

# Add colourbar to image
fig.add_colorbar()
fig.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot C Band contours
sigma=1.1
fig.show_contour(image_mJy, levels=[-3*sigma,
    25*sigma, 50*sigma, 75*sigma,
    100*sigma,125*sigma, 150*sigma, 175*sigma],
	linewidths=0.5, colors='white', overlap=True)

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
fig.beam.set_color(beamcolor)
fig.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
fig.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig.axis_labels.set_font(size=labelfontsize)
fig.tick_labels.set_font(size=labelfontsize)
fig.colorbar.set_axis_label_font(size=labelfontsize)
fig.colorbar.set_font(size=labelfontsize)
fig.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig.frame.set_color(framecolor)
#Set tick colour to black
fig.ticks.set_color(tickcolor)

# Add title and labels for components
fig.add_label(0.33, 0.28, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig.add_label(0.67, 0.71, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig.add_label(0.5, 0.92, 'Band 7 (336 GHz)', relative=True, color=labelcolor, size=labelfontsize)

fig.save(final_image, dpi=500)

