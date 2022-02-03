from radmc3dPy.image import *
from matplotlib import cm
#a=readImage("K_Band_grid_1000_zmin_8au_w0_1au_v0_100.out")
a=readImage("image.out")
#a = a.imConv(dpc=140.,psfType='gauss',fwhm=[0.081,0.034],pa=-28.29)
a = a.imConv(dpc=140.,psfType='gauss',fwhm=[0.104,0.081],pa=-23.26)
plotImage(a,log=False,cmap=cm.hot,bunit='jy/beam',dpc=140,au=False,arcsec=True)
