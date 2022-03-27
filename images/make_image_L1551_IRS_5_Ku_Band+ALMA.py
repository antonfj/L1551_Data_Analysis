import numpy as np
import matplotlib
import os
import aplpy
from astropy.io import fits

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
alma_image = 'L1551_IRS_5_ALMA_Band_4_shifted_to_VLA_epoch.fits'
alma_image_mJy = 'L1551_IRS_5_ALMA_Band_4_shifted_to_VLA_epoch.mJy.fits'
ku_band_image = 'L1551_IRS_5_Ku_Band_VLA_r_0_5_sc1.cropped.fits'
ku_band_image_mJy = 'L1551_IRS_5_Ku_Band_VLA_r_0_5_sc1.cropped.mJy.fits'

final_image = 'L1551_IRS_5_ALMA+Ku_Band_VLA.pdf'

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

convert_to_mJy(alma_image,alma_image_mJy)
convert_to_mJy(ku_band_image,ku_band_image_mJy)

titlefontsize=12
labelfontsize=8
cmap='jet'
framecolor='white'
tickcolor='white'
labelcolor='white'
beamcolor='white'

fig = aplpy.FITSFigure(alma_image_mJy,figsize=[6.,5.])
fig.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig.show_colorscale(stretch='linear')
fig.show_colorscale(vmin=-0.2, vmax=13, cmap='jet')

# Add colourbar to image
fig.add_colorbar()
fig.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot Ku Band contours
sigma=4.3e-3
fig.show_contour(ku_band_image_mJy, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	6*sigma, 7*sigma, 8*sigma, 9*sigma, 10*sigma, 20*sigma, 40*sigma, 60*sigma,
	80*sigma, 100*sigma, 120*sigma, 140*sigma, 160*sigma, 180*sigma], colors='white', linewidths=0.5, overlap=True)

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
fig.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig.add_beam()
# Set to axes and PA of beam of 15 GHz image
fig.beam.set_major(4.44e-5)    # 0.16 arcsecs
fig.beam.set_minor(3.33e-5)    # 0.12 arcsecs
fig.beam.set_angle(-60)
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

#Set frame colour to black
fig.frame.set_color(framecolor)
#Set tick colour to black
fig.ticks.set_color(tickcolor)

fig.save(final_image, dpi=500)

