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
c_band_image = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.fits'
c_band_image_mJy = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.mJy.fits'
ku_band_image = 'L1551_IRS_5_Ku_Band_VLA_r_0_5_sc1.cropped.fits'
ku_band_image_mJy = 'L1551_IRS_5_Ku_Band_VLA_r_0_5_sc1.cropped.mJy.fits'

final_image = 'L1551_IRS_5_C_Band_Ku_Band_combined.pdf'

fig = matplotlib.pyplot.figure(figsize=[7.2,2.8])

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
convert_to_mJy(c_band_image,c_band_image_mJy)

titlefontsize=12
labelfontsize=8
cmap='jet'
framecolor='white'
tickcolor='white'
labelcolor='white'
beamcolor='white'


#########################################################################################
########################### C Band image ################################################
#########################################################################################

fig1 = aplpy.FITSFigure(c_band_image_mJy, figure=fig, subplot=[0.14, 0.08, 0.36, 0.92])
fig1.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig1.show_colorscale(stretch='linear')
fig1.show_colorscale(vmin=-0.03, vmax=0.3, cmap='jet')

# Add colourbar to image
fig1.add_colorbar()
fig1.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot C Band contours
sigma=11e-3
fig1.show_contour(c_band_image_mJy, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	6*sigma, 7*sigma, 8*sigma, 9*sigma, 10*sigma, 15*sigma, 20*sigma, 25*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Coordinates of sources
north_coord = np.array([67.8923542,18.1346203])
south_coord = np.array([67.8923708,18.1345203])

# Jet axis of sources
jet_axis_north = 67.0 * (np.pi/180.0)
jet_axis_south = 55.0 * (np.pi/180.0)

# Plot line along jet axes of both sources
ra_length = 2.0 / (60*60)	# Set length of line along ra axis
dec_length_north = ra_length / np.tan(jet_axis_north)
dec_length_south = ra_length / np.tan(jet_axis_south)

jet_coords = np.array([[ [north_coord[0] - ra_length/2, north_coord[0] + ra_length/2],
	[north_coord[1] - dec_length_north/2, north_coord[1] + dec_length_north/2] ],
	[ [south_coord[0] - ra_length/2, south_coord[0] + ra_length/2],
	[south_coord[1] - dec_length_south/2, south_coord[1] + dec_length_south/2] ] ])
fig1.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig1.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig1.add_beam()
# Set to axes and PA of beam of 5 GHz image
fig1.beam.set_major(6.11e-5)    # 0.22 arcsecs
fig1.beam.set_minor(3.06e-5)    # 0.11 arcsecs
fig1.beam.set_angle(-44)
fig1.beam.set_color(beamcolor)
fig1.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
fig1.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig1.axis_labels.set_font(size=labelfontsize)
fig1.tick_labels.set_font(size=labelfontsize)
fig1.colorbar.set_axis_label_font(size=labelfontsize)
fig1.colorbar.set_font(size=labelfontsize)
fig1.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig1.frame.set_color(framecolor)
#Set tick colour to black
fig1.ticks.set_color(tickcolor)


#########################################################################################
########################### Ku Band image ###############################################
#########################################################################################

fig2 = aplpy.FITSFigure(ku_band_image_mJy, figure=fig, subplot=[0.58, 0.08, 0.36, 0.92])
fig2.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig2.show_colorscale(stretch='linear')
fig2.show_colorscale(vmin=-0.01, vmax=0.7, cmap='jet')

# Add colourbar to image
fig2.add_colorbar()
fig2.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot Ku Band contours
sigma=4.3e-3
fig2.show_contour(ku_band_image_mJy, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	6*sigma, 7*sigma, 8*sigma, 9*sigma, 10*sigma, 20*sigma, 40*sigma, 60*sigma,
	80*sigma, 100*sigma, 120*sigma, 140*sigma, 160*sigma, 180*sigma], colors='white', linewidths=0.5, overlap=True)

# Plot line along jet axes of both sources
fig2.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig2.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig2.add_beam()
# Set to axes and PA of beam of 15 GHz image
fig2.beam.set_major(4.44e-5)    # 0.16 arcsecs
fig2.beam.set_minor(3.33e-5)    # 0.12 arcsecs
fig2.beam.set_angle(-60)
fig2.beam.set_color(beamcolor)
fig2.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
fig2.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig2.axis_labels.set_font(size=labelfontsize)
fig2.tick_labels.set_font(size=labelfontsize)
fig2.colorbar.set_axis_label_font(size=labelfontsize)
fig2.colorbar.set_font(size=labelfontsize)
fig2.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig2.frame.set_color(framecolor)
#Set tick colour to black
fig2.ticks.set_color(tickcolor)

# Remove RA labels
fig2.axis_labels.hide_y()
fig2.tick_labels.hide_y()

fig.savefig(final_image, dpi=500)