# test space for learning to use pyds9

import pyds9
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
import numpy as np

# check if ds9 is accesible
if not pyds9.ds9_targets():
    input('DS9 target not found. Please start/restart DS9, then press enter')



# show targets
print('ds9 target instance')
print(pyds9.ds9_targets())

# pull the region data
region_data = get_regions_from_ds9()
region1 = region_data[0]

norm = ImageNormalize(stretch=SqrtStretch())
plt.figure()
plt.imshow(region1.data, norm=norm, origin='lower', cmap='viridis')
plt.show()


# need to have an instance of ds9 running

# define the ds9 object.  Calls class DS9
ds9 = pyds9.DS9()

# make a new, seperate instance of ds9
ds9display = pyds9.DS9(target='display', start='-title display')

# open file
filename = '/home/lee/Documents/k4m_160319_101212_ori.fits.fz[im2]'
ds9.set('file ' + filename)


# show the open file in the target instance of ds9
print(ds9.get('file'))


# show the xpa id of target instance
print(ds9.access())

# show the region(s)
print(ds9.get('regions getinfo'))

# set region format, in this case to image
ds9.set('regions system image')

# show only the regions that have been selected
ds9.get('regions selected')

# select all defined regions
ds9.set('regions select all')

# open fits file as mosaic
ds9.set('mosaicimage ' + filename)

# open fits file as multiple frames
ds9.set('multiframe ' + filename)

# accessing data slices directly
data_str = ds9.get('data image 940 1092 8 6 no')

# convert to 1d array
# data_1d = np.fromstring(data_str, sep='\n')
# above does not work

# convert to original shape

