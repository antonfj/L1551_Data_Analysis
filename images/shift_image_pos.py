from astropy.io import fits

image_file='L1551_IRS5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut.cropped.fits'
output_file='L1551_IRS5_C_Band_E-MERLIN_VLA_combined_r_0.5_uv_cut_shifted_to_VLA_pos.cropped.fits'

image_data, image_header = fits.getdata(image_file, header=True)
#print(image_header)

# RA and Dec shifts in milliarcsec
ra_shift = -34.7
dec_shift = 92.2

# Convert to degrees
ra_shift = ra_shift/(1000*60*60)
dec_shift = dec_shift/(1000*60*60)

print("RA Shift (deg): ", ra_shift)
print("Dec Shift (deg): ", dec_shift)

print("RA value at centre of image: ", image_header['CRVAL1'])
print("Dec value at centre of image: ", image_header['CRVAL2'])

# Shift positions of RA and Dec values in image
image_header['CRVAL1'] -= ra_shift
image_header['CRVAL2'] -= dec_shift

print("New RA value at centre of image: ", image_header['CRVAL1'])
print("New Dec value at centre of image: ", image_header['CRVAL2'])

# Write image to new FITS file
fits.writeto(output_file, image_data, image_header)
