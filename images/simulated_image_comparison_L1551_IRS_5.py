import numpy as np
import matplotlib
import os
import aplpy
from astropy.io import fits

# Set tick directions manually
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Images to use in the final image
simulated_c_image = 'L1551_IRS_5_C_Band_5_GHz_simulated.fits'
simulated_c_image_mJy = 'L1551_IRS_5_C_Band_5_GHz_simulated.mJy.fits'
observed_c_image = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.fits'
observed_c_image_mJy = 'L1551_IRS_5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.mJy.fits'

simulated_x_image = 'L1551_IRS_5_X_Band_simulated.fits'
simulated_x_image_mJy = 'L1551_IRS_5_X_Band_simulated.mJy.fits'
observed_x_image = 'L1551_IRS5_X_Band_VLA_initial_r_+0.5.fits'
observed_x_image_mJy = 'L1551_IRS5_X_Band_VLA_initial_r_+0.5.mJy.fits'

final_image = 'L1551_IRS_5_simulated_image_comparison.pdf'

# Convert images to mJy
# Check if mJy files exist

def convert_to_mJy(image, image_mJy):
    if not os.path.isfile(image_mJy):
        # Copy to new file
        os.system('cp ' + image + ' ' + image_mJy)

        # Opens FITS image
        hdulist = fits.open(image_mJy, mode='update')
        fitsdata = hdulist[0].data[0,:]
        
        # Converts to uJy
        fitsdata = fitsdata*1.e3

        hdulist[0].data[0,:] = fitsdata
        hdulist.flush()

convert_to_mJy(simulated_c_image, simulated_c_image_mJy)
convert_to_mJy(observed_c_image, observed_c_image_mJy)

convert_to_mJy(simulated_x_image, simulated_x_image_mJy)
convert_to_mJy(observed_x_image, observed_x_image_mJy)

fig = matplotlib.pyplot.figure(figsize=[7,5.])

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
############### Observed C Band image #################### 
###################################################
fig1 = aplpy.FITSFigure(observed_c_image_mJy, figure=fig, subplot=[0.12, 0.53, 0.4, 0.42])
fig1.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)

fig1.show_colorscale(stretch='linear')
fig1.show_colorscale(vmin=-2e-2, vmax=3.2e-1, cmap=cmap)

# Add colourbar to image
fig1.add_colorbar()
fig1.colorbar.set_font(size=labelfontsize)
fig1.colorbar.set_axis_label_font(size=labelfontsize)
fig1.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Set title
fig1.set_title('5 GHz Observation', size=titlefontsize)

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

# Remove RA/Dec labels
fig1.axis_labels.hide_x()
fig1.tick_labels.hide_x()

###################################################
############### Simulated C Band image ################### 
###################################################
fig2 = aplpy.FITSFigure(simulated_c_image_mJy, figure=fig, subplot=[0.12, 0.06, 0.4, 0.42])
fig2.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)
fig2.show_colorscale(stretch='linear')
fig2.show_colorscale(vmin=-2e-2, vmax=3.2e-1, cmap=cmap)

# Add colourbar to image
fig2.add_colorbar()
fig2.colorbar.set_font(size=labelfontsize)
fig2.colorbar.set_axis_label_font(size=labelfontsize)
fig2.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Set title
fig2.set_title('5 GHz Model', size=titlefontsize)

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

###################################################
############### Observed X Band image #################### 
###################################################
fig3 = aplpy.FITSFigure(observed_x_image_mJy, figure=fig, subplot=[0.55, 0.53, 0.4, 0.42])
fig3.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)

fig3.show_colorscale(stretch='linear')
fig3.show_colorscale(vmin=-2e-2, vmax=7e-1, cmap=cmap)

# Add colourbar to image
fig3.add_colorbar()
fig3.colorbar.set_font(size=labelfontsize)
fig3.colorbar.set_axis_label_font(size=labelfontsize)
fig3.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Set title
fig3.set_title('10 GHz Observation', size=titlefontsize)

# Adds synthesis beam in bottom left corner
#fig.add_beam(major=0.0055, minor=0.0055, angle=90.0)
#fig1.add_beam()
#fig1.beam.set_color(beamcolor)
#fig1.beam.set_pad(1.0)

# Plot marker at position of both sources
fig3.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Set font size of labels
fig3.axis_labels.set_font(size=labelfontsize)
fig3.tick_labels.set_font(size=labelfontsize)

#Set frame colour
fig3.frame.set_color(framecolor)
#Set tick colour
fig3.ticks.set_color(tickcolor)

# Remove RA/Dec labels
fig3.axis_labels.hide_x()
fig3.tick_labels.hide_x()
fig3.axis_labels.hide_y()
fig3.tick_labels.hide_y()

###################################################
############### Simulated C Band image ################### 
###################################################
fig4 = aplpy.FITSFigure(simulated_x_image_mJy, figure=fig, subplot=[0.55, 0.06, 0.4, 0.42])
fig4.recenter(67.8923708,18.1345203, width=0.2e-3, height=0.2e-3)
fig4.show_colorscale(stretch='linear')
fig4.show_colorscale(vmin=-2e-2, vmax=7e-1, cmap=cmap)

# Add colourbar to image
fig4.add_colorbar()
fig4.colorbar.set_font(size=labelfontsize)
fig4.colorbar.set_axis_label_font(size=labelfontsize)
fig4.colorbar.set_axis_label_text(r'Flux (mJy/beam)')

# Set title
fig4.set_title('10 GHz Model', size=titlefontsize)

# Adds synthesis beam in bottom left corner
#fig.add_beam(major=0.0055, minor=0.0055, angle=90.0)
#fig4.add_beam()
#fig4.beam.set_color(beamcolor)
#fig4.beam.set_pad(1.0)

# Plot marker at position of both sources
fig4.show_markers([north_coord[0], south_coord[0]], [north_coord[1], south_coord[1]], marker='+', 
                 facecolor='black', edgecolor='black', zorder=10)

# Set font size of labels
fig4.axis_labels.set_font(size=labelfontsize)
fig4.tick_labels.set_font(size=labelfontsize)

# Remove Dec labels
#fig4.axis_labels.hide_y()
#fig4.tick_labels.hide_y()

#Set frame colour
fig4.frame.set_color(framecolor)
#Set tick colour
fig4.ticks.set_color(tickcolor)

# Remove RA/Dec labels
fig4.axis_labels.hide_y()
fig4.tick_labels.hide_y()

fig.savefig(final_image, dpi=500)

matplotlib.pyplot.show()