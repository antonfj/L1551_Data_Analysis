import numpy as np
import matplotlib
import os
import aplpy
from astropy.io import fits

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
simulated_image = 'L1551_IRS_5_C_Band_5_GHz_simulated.fits'
simulated_image_mJy = 'L1551_IRS_5_C_Band_5_GHz_simulated.mJy.fits'
observed_image = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.fits'
observed_image_mJy = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.mJy.fits'

final_image = 'L1551_IRS_5_C_Band_simulated_image_comparison.pdf'

# Convert images to mJy
# Check if mJy files exist
if not os.path.isfile(simulated_image_mJy):
    # Copy to new file
    os.system('cp ' + simulated_image + ' ' + simulated_image_mJy)
    
    # Opens FITS image
    hdulist = fits.open(simulated_image_mJy, mode='update')
    fitsdata = hdulist[0].data

    # Converts to uJy
    fitsdata = fitsdata*1.e3

    hdulist[0].data = fitsdata
    hdulist.flush()

if not os.path.isfile(observed_image_mJy):
    # Copy to new file
    os.system('cp ' + observed_image + ' ' + observed_image_mJy)
    
    # Opens FITS image
    hdulist = fits.open(observed_image_mJy, mode='update')
    fitsdata = hdulist[0].data[:,0]

    # Converts to uJy
    fitsdata = fitsdata*1.e3

    hdulist[0].data[:,0] = fitsdata
    hdulist.flush()


fig = matplotlib.pyplot.figure(figsize=[3.5,5.])

titlefontsize=8
labelfontsize=6
cmap='jet'
framecolor='black'
tickcolor='white'
labelcolor='white'
beamcolor='white'

# Coordinates of sources
north_coord = np.array([67.8923542,18.1346203])
south_coord = np.array([67.8923708,18.1345203])

###################################################
############### Observed image #################### 
###################################################
fig1 = aplpy.FITSFigure(observed_image_mJy, figure=fig, subplot=[0.13, 0.56, 0.85, 0.4])
fig1.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)

fig1.show_colorscale(stretch='linear')
fig1.show_colorscale(vmin=-2e-2, vmax=3.2e-1, cmap=cmap)

# Add colourbar to image
fig1.add_colorbar()
fig1.colorbar.set_font(size=labelfontsize)
fig1.colorbar.set_axis_label_font(size=labelfontsize)
fig1.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Set title
fig1.set_title('Observation', size=titlefontsize)

# Adds synthesis beam in bottom left corner
#fig.add_beam(major=0.0055, minor=0.0055, angle=90.0)
#fig1.add_beam()
#fig1.beam.set_color(beamcolor)
#fig1.beam.set_pad(1.0)

# Plot marker at position of both sources
fig1.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Set font size of labels
fig1.axis_labels.set_font(size=labelfontsize)
fig1.tick_labels.set_font(size=labelfontsize)

#Set frame colour
fig1.frame.set_color(framecolor)
#Set tick colour
fig1.ticks.set_color(tickcolor)

###################################################
############### Simulated image ################### 
###################################################
fig2 = aplpy.FITSFigure(simulated_image_mJy, figure=fig, subplot=[0.13, 0.06, 0.85, 0.4])
fig2.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)
fig2.show_colorscale(stretch='linear')
fig2.show_colorscale(vmin=-2e-2, vmax=3.2e-1, cmap=cmap)

# Add colourbar to image
fig2.add_colorbar()
fig2.colorbar.set_font(size=labelfontsize)
fig2.colorbar.set_axis_label_font(size=labelfontsize)
fig2.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Set title
fig2.set_title('Model', size=titlefontsize)

# Adds synthesis beam in bottom left corner
#fig.add_beam(major=0.0055, minor=0.0055, angle=90.0)
#fig2.add_beam()
#fig2.beam.set_color(beamcolor)
#fig2.beam.set_pad(1.0)

# Plot marker at position of both sources
fig2.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Set font size of labels
fig2.axis_labels.set_font(size=labelfontsize)
fig2.tick_labels.set_font(size=labelfontsize)

# Remove Dec labels
#fig2.axis_labels.hide_y()
#fig2.tick_labels.hide_y()

#Set frame colour
fig2.frame.set_color(framecolor)
#Set tick colour
fig2.ticks.set_color(tickcolor)

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