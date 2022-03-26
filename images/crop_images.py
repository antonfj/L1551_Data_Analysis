import os

full_image='L1551_IRS5_X_Band_VLA_widefield_initial_r_+0.5.pbcor.fits'
cropped_image='L1551_NE_X_Band_VLA_cropped_r_+0.5.pbcor.image'
cropped_image_fits='L1551_NE_X_Band_VLA_cropped_r_+0.5.pbcor.fits'

os.system('rm ' + cropped_image_fits)

imsubimage(imagename=full_image, outfile=cropped_image, box="850, 4330, 1250, 4730")

exportfits(imagename=cropped_image, fitsimage=cropped_image_fits)

rmtables(cropped_image)
