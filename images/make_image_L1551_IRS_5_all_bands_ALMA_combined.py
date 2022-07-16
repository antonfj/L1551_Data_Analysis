import numpy as np
import matplotlib
import os
import aplpy
from astropy.io import fits

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
band_3_image_mJy = 'L1551_IRS_5_ALMA_Band_3_cropped_shifted_to_VLA_pos.image.pbcor.mJy.fits'
band_4_image_mJy = 'L1551_IRS_5_ALMA_Band_4_shifted_to_VLA_epoch.mJy.fits'
band_6_image_mJy = 'L1551_IRS_5_ALMA_Band_6_shifted.mJy.fits'
band_7_image_mJy = 'L1551_IRS_5_ALMA_Band_7_shifted_to_VLA_pos.cont.I.pbcor.mJy.fits'

final_image = 'L1551_IRS_5_all_bands_ALMA_combined.pdf'

fig = matplotlib.pyplot.figure(figsize=[7, 5.2])

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

titlefontsize=12
labelfontsize=8
cmap='jet'
framecolor='white'
tickcolor='white'
labelcolor='red'
beamcolor='white'


#########################################################################################
########################### Band 3 image ################################################
#########################################################################################

fig1 = aplpy.FITSFigure(band_3_image_mJy, figure=fig, subplot=[0.145, 0.52, 0.36, 0.5])
fig1.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig1.show_colorscale(stretch='linear')
fig1.show_colorscale(vmin=-150e-3, vmax=15, cmap='jet')

# Add colourbar to image
fig1.add_colorbar()
fig1.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot contours
sigma=50e-3
fig1.show_contour(band_3_image_mJy, levels=[-3*sigma, 50*sigma,
    100*sigma, 150*sigma, 200*sigma, 250*sigma, 300*sigma],
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
fig1.beam.set_color(beamcolor)
fig1.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
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

# Add title and labels for components
fig1.add_label(0.42, 0.33, 'S', relative=True, color=labelcolor, size=labelfontsize)
fig1.add_label(0.59, 0.65, 'N', relative=True, color=labelcolor, size=labelfontsize)
fig1.add_label(0.5, 0.92, '(a) Band 3 (93 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
fig1.axis_labels.hide_x()
fig1.tick_labels.hide_x()

#########################################################################################
########################### Band 4 image ###############################################
#########################################################################################

fig2 = aplpy.FITSFigure(band_4_image_mJy, figure=fig, subplot=[0.575, 0.52, 0.36, 0.5])
fig2.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig2.show_colorscale(stretch='linear')
fig2.show_colorscale(vmin=-250e-3, vmax=13, cmap='jet')

# Add colourbar to image
fig2.add_colorbar()
fig2.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot contours
sigma=76e-3
fig2.show_contour(band_4_image_mJy, levels=[-3*sigma,
	25*sigma, 50*sigma, 75*sigma, 100*sigma,125*sigma, 150*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig2.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig2.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig2.add_beam()
fig2.beam.set_color(beamcolor)
fig2.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
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

# Add title and labels for components
fig2.add_label(0.44, 0.32, 'S', relative=True, color=labelcolor, size=labelfontsize)
fig2.add_label(0.59, 0.65, 'N', relative=True, color=labelcolor, size=labelfontsize)
fig2.add_label(0.5, 0.92, '(b) Band 4 (153 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
fig2.axis_labels.hide_x()
fig2.tick_labels.hide_x()
fig2.axis_labels.hide_y()
fig2.tick_labels.hide_y()

#########################################################################################
########################### Band 6 image ###############################################
#########################################################################################

fig3 = aplpy.FITSFigure(band_6_image_mJy, figure=fig, subplot=[0.145, 0.05, 0.36, 0.5])
fig3.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig3.show_colorscale(stretch='linear')
fig3.show_colorscale(vmin=-3, vmax=130, cmap='jet')

# Add colourbar to image
fig3.add_colorbar()
fig3.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot contours
sigma=0.75
fig3.show_contour(band_6_image_mJy, levels=[-3*sigma,
    25*sigma, 50*sigma, 75*sigma,
    100*sigma,125*sigma, 150*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig3.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig3.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig3.add_beam()
fig3.beam.set_color(beamcolor)
fig3.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
fig3.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig3.axis_labels.set_font(size=labelfontsize)
fig3.tick_labels.set_font(size=labelfontsize)
fig3.colorbar.set_axis_label_font(size=labelfontsize)
fig3.colorbar.set_font(size=labelfontsize)
fig3.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig3.frame.set_color(framecolor)
#Set tick colour to black
fig3.ticks.set_color(tickcolor)

# Add title and labels for components
fig3.add_label(0.44, 0.32, 'S', relative=True, color=labelcolor, size=labelfontsize)
fig3.add_label(0.59, 0.65, 'N', relative=True, color=labelcolor, size=labelfontsize)
fig3.add_label(0.5, 0.92, '(c) Band 6 (225 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
#fig3.axis_labels.hide_x()
#fig3.tick_labels.hide_x()
#fig3.axis_labels.hide_y()
#fig3.tick_labels.hide_y()

#########################################################################################
########################### Band 7 image ###############################################
#########################################################################################

fig4 = aplpy.FITSFigure(band_7_image_mJy, figure=fig, subplot=[0.575, 0.05, 0.36, 0.5])
fig4.recenter(67.8923625, 18.1345753, width=0.7e-3, height=0.7e-3)

fig4.show_colorscale(stretch='linear')
fig4.show_colorscale(vmin=-3, vmax=150, cmap='jet')

# Add colourbar to image
fig4.add_colorbar()
fig4.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot contours
sigma=1.1
fig4.show_contour(band_7_image_mJy, levels=[-3*sigma,
    25*sigma, 50*sigma, 75*sigma,
    100*sigma,125*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig4.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig4.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig4.add_beam()
fig4.beam.set_color(beamcolor)
fig4.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
fig4.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig4.axis_labels.set_font(size=labelfontsize)
fig4.tick_labels.set_font(size=labelfontsize)
fig4.colorbar.set_axis_label_font(size=labelfontsize)
fig4.colorbar.set_font(size=labelfontsize)
fig4.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig4.frame.set_color(framecolor)
#Set tick colour to black
fig4.ticks.set_color(tickcolor)

# Add title and labels for components
fig4.add_label(0.44, 0.32, 'S', relative=True, color=labelcolor, size=labelfontsize)
fig4.add_label(0.59, 0.65, 'N', relative=True, color=labelcolor, size=labelfontsize)
fig4.add_label(0.5, 0.92, '(d) Band 7 (336 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
#fig4.axis_labels.hide_x()
#fig4.tick_labels.hide_x()
fig4.axis_labels.hide_y()
fig4.tick_labels.hide_y()

###############################################################################
######################## Make final combined image ############################
###############################################################################

fig.savefig(final_image, dpi=500)