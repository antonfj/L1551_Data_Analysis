import numpy as np
import matplotlib
import os
import aplpy
from astropy.io import fits

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
c_band_image_mJy = 'L1551_NE_C_Band_VLA_initial_r_+0.5.pbcor.cropped.mJy.fits'
x_band_image_mJy = 'L1551_NE_X_Band_VLA_cropped_r_+0.5.pbcor.mJy.fits'
ku_band_lower_image_mJy = 'L1551_NE_Ku_Band_VLA_12-15_GHz_r_+0.5.mJy.fits'
ku_band_upper_image_mJy = 'L1551_NE_Ku_Band_VLA_15-18_GHz_r_+0.5.mJy.fits'
k_band_lower_image_mJy = 'L1551_NE_K_Band_VLA_18-22_GHz_r_+0.5.mJy.fits'
k_band_upper_image_mJy = 'L1551_NE_K_Band_VLA_22-26_GHz_r_+0.5.mJy.fits'

final_image = 'L1551_NE_all_bands_C_to_K_combined.pdf'

fig = matplotlib.pyplot.figure(figsize=[7, 7.5])

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
########################### C Band image ################################################
#########################################################################################

fig1 = aplpy.FITSFigure(c_band_image_mJy, figure=fig, subplot=[0.15, 0.68, 0.36, 0.3])
fig1.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig1.show_colorscale(stretch='linear')
fig1.show_colorscale(vmin=-5e-2, vmax=1.6e-1, cmap='jet')

# Add colourbar to image
fig1.add_colorbar()
fig1.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot C Band contours
sigma=16e-3
fig1.show_contour(c_band_image_mJy, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	6*sigma, 7*sigma, 8*sigma, 9*sigma],
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
fig1.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'])

# Plot marker at position of both sources
fig1.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig1.add_beam()
# Set to axes and PA of beam of 5 GHz image
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

# Add title and labels for components
fig1.add_label(0.38, 0.33, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig1.add_label(0.62, 0.66, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig1.add_label(0.5, 0.92, '(a) C Band (5 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
fig1.axis_labels.hide_x()
fig1.tick_labels.hide_x()

#########################################################################################
########################### X Band image ###############################################
#########################################################################################

fig2 = aplpy.FITSFigure(x_band_image_mJy, figure=fig, subplot=[0.58, 0.68, 0.36, 0.3])
fig2.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig2.show_colorscale(stretch='linear')
fig2.show_colorscale(vmin=-8e-2, vmax=3e-1, cmap='jet')

# Add colourbar to image
fig2.add_colorbar()
fig2.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot X Band contours
sigma=26e-3
fig2.show_contour(x_band_image_mJy, levels=[-3*sigma, 3*sigma, 4*sigma, 5*sigma,
	6*sigma, 7*sigma, 8*sigma, 9*sigma, 10*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig2.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig2.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig2.add_beam()
# Set to axes and PA of beam of 15 GHz image
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

# Add title and labels for components
fig2.add_label(0.38, 0.34, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig2.add_label(0.62, 0.65, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig2.add_label(0.5, 0.92, '(b) X Band (10 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
fig2.axis_labels.hide_x()
fig2.tick_labels.hide_x()
fig2.axis_labels.hide_y()
fig2.tick_labels.hide_y()

#########################################################################################
########################### Ku Band lower image ###############################################
#########################################################################################

fig3 = aplpy.FITSFigure(ku_band_lower_image_mJy, figure=fig, subplot=[0.15, 0.365, 0.36, 0.3])
fig3.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig3.show_colorscale(stretch='linear')
fig3.show_colorscale(vmin=-2e-2, vmax=3.5e-1, cmap='jet')

# Add colourbar to image
fig3.add_colorbar()
fig3.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot Ku Band contours
sigma=6e-3
fig3.show_contour(ku_band_lower_image_mJy, levels=[-3*sigma, 3*sigma,
	10*sigma, 20*sigma, 30*sigma, 40*sigma, 50*sigma, 60*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig3.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig3.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig3.add_beam()
# Set to axes and PA of beam of 15 GHz image
fig3.beam.set_color(beamcolor)
fig3.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
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
fig3.add_label(0.38, 0.34, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig3.add_label(0.62, 0.65, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig3.add_label(0.5, 0.92, '(c) Ku Band (13.5 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
fig3.axis_labels.hide_x()
fig3.tick_labels.hide_x()
#fig3.axis_labels.hide_y()
#fig3.tick_labels.hide_y()

#########################################################################################
########################### Ku Band upper image ###############################################
#########################################################################################

fig4 = aplpy.FITSFigure(ku_band_upper_image_mJy, figure=fig, subplot=[0.58, 0.365, 0.36, 0.3])
fig4.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig4.show_colorscale(stretch='linear')
fig4.show_colorscale(vmin=-2e-2, vmax=4.5e-1, cmap='jet')

# Add colourbar to image
fig4.add_colorbar()
fig4.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot Ku Band contours
sigma=6e-3
fig4.show_contour(ku_band_upper_image_mJy, levels=[-3*sigma, 3*sigma,
	10*sigma, 20*sigma, 30*sigma, 40*sigma, 50*sigma, 60*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig4.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig4.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig4.add_beam()
fig4.beam.set_color(beamcolor)
fig4.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
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
fig4.add_label(0.38, 0.34, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig4.add_label(0.62, 0.65, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig4.add_label(0.5, 0.92, '(d) Ku Band (17.5 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA/Dec labels
fig4.axis_labels.hide_x()
fig4.tick_labels.hide_x()
fig4.axis_labels.hide_y()
fig4.tick_labels.hide_y()

#########################################################################################
########################### K Band lower image ###############################################
#########################################################################################

fig5 = aplpy.FITSFigure(k_band_lower_image_mJy, figure=fig, subplot=[0.15, 0.05, 0.36, 0.3])
fig5.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig5.show_colorscale(stretch='linear')
fig5.show_colorscale(vmin=-3e-2, vmax=4.5e-1, cmap='jet')

# Add colourbar to image
fig5.add_colorbar()
fig5.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot Ku Band contours
sigma=11e-3
fig5.show_contour(k_band_lower_image_mJy, levels=[-3*sigma, 3*sigma, 5*sigma,
	10*sigma, 15*sigma, 20*sigma, 25*sigma, 30*sigma, 35*sigma, 40*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig5.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig5.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig5.add_beam()
# Set to axes and PA of beam of 15 GHz image
fig5.beam.set_major(4.44e-5)    # 0.16 arcsecs
fig5.beam.set_minor(3.33e-5)    # 0.12 arcsecs
fig5.beam.set_angle(-60)
fig5.beam.set_color(beamcolor)
fig5.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
fig5.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig5.axis_labels.set_font(size=labelfontsize)
fig5.tick_labels.set_font(size=labelfontsize)
fig5.colorbar.set_axis_label_font(size=labelfontsize)
fig5.colorbar.set_font(size=labelfontsize)
fig5.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig5.frame.set_color(framecolor)
#Set tick colour to black
fig5.ticks.set_color(tickcolor)

# Add title and labels for components
fig5.add_label(0.38, 0.34, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig5.add_label(0.62, 0.65, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig5.add_label(0.5, 0.92, '(e) K Band (20 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA labels
#fig3.axis_labels.hide_y()
#fig3.tick_labels.hide_y()

#########################################################################################
########################### K Band upper image ###############################################
#########################################################################################

fig6 = aplpy.FITSFigure(k_band_upper_image_mJy, figure=fig, subplot=[0.58, 0.05, 0.36, 0.3])
fig6.recenter(67.9354, 18.14205, width=0.7e-3, height=0.7e-3)

fig6.show_colorscale(stretch='linear')
fig6.show_colorscale(vmin=-3e-2, vmax=4e-1, cmap='jet')

# Add colourbar to image
fig6.add_colorbar()
fig6.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Overplot Ku Band contours
sigma=11e-3
fig6.show_contour(k_band_upper_image_mJy, levels=[-3*sigma, 3*sigma, 5*sigma,
	10*sigma, 15*sigma, 20*sigma, 25*sigma, 30*sigma, 35*sigma],
	linewidths=0.5, colors='white', overlap=True)

# Plot line along jet axes of both sources
fig6.show_lines(jet_coords, color=['white', 'white'], linestyle=['dashed', 'dashed'], linewidths=1)

# Plot marker at position of both sources
fig6.show_markers([A_coord[0], B_coord[0]], [A_coord[1], B_coord[1]], marker='+', facecolor='black', edgecolor='black', zorder=10)

# Adds synthesis beam in bottom left corner
fig6.add_beam()
# Set to axes and PA of beam of 15 GHz image
fig6.beam.set_major(4.44e-5)    # 0.16 arcsecs
fig6.beam.set_minor(3.33e-5)    # 0.12 arcsecs
fig6.beam.set_angle(-60)
fig6.beam.set_color(beamcolor)
fig6.beam.set_pad(1.0)

# Add scale bar
length = 20     # length in au
angular_length = (length / (206265 * 147)) * (360/(2*np.pi))
print(angular_length)
fig6.add_scalebar(angular_length, label='20 AU', corner='bottom right', color=beamcolor)

# Set font size of labels
fig6.axis_labels.set_font(size=labelfontsize)
fig6.tick_labels.set_font(size=labelfontsize)
fig6.colorbar.set_axis_label_font(size=labelfontsize)
fig6.colorbar.set_font(size=labelfontsize)
fig6.scalebar.set_font_size(size=labelfontsize)

#Set frame colour to black
fig6.frame.set_color(framecolor)
#Set tick colour to black
fig6.ticks.set_color(tickcolor)

# Add title and labels for components
fig6.add_label(0.38, 0.34, 'A', relative=True, color=labelcolor, size=labelfontsize)
fig6.add_label(0.62, 0.65, 'B', relative=True, color=labelcolor, size=labelfontsize)
fig6.add_label(0.5, 0.92, '(f) K Band (24 GHz)', relative=True, color=labelcolor, size=labelfontsize)

# Remove RA labels
fig6.axis_labels.hide_y()
fig6.tick_labels.hide_y()

fig.savefig(final_image, dpi=500)