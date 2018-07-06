# test space for displaying regions from SAOImage DS9

import pyds9 as ds9
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt

# show targets
print('ds9 target instance')
print(ds9.ds9_targets())

# pull the region data
from get_regions import get_regions

region_data = get_regions()
region1 = region_data[0]

norm = ImageNormalize(stretch=SqrtStretch())
plt.figure()
plt.imshow(region1.data, norm=norm, origin='lower', cmap='viridis')
plt.show()

# call the fits data externally
from astropy.io import fits

# open fits file, best practice
file_name = '/home/lee/Documents/k4m_160319_101212_ori.fits.fz'
with fits.open(file_name) as hdu:
    hdu.info()
    data_im1 = hdu[1].data

# slice the data
ymin = region1.ymin
ymax = region1.ymax
xmin = region1.xmin
xmax = region1.xmax
Object1_Data = data_im1[ymin:ymax,xmin:xmax]

# displace the result
norm = ImageNormalize(stretch=SqrtStretch())
plt.figure()
plt.imshow(region1.data, norm=norm, origin='lower', cmap='viridis')
plt.show()
