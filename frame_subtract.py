"""
A function for subtracting (and adding) frames from fits files
"""
import numpy as np
import pyds9
from astropy.io import fits

import timing
def frame_subtract(minuend, subtrahend):
    difference = [minuend.data.astype('float64', casting='safe') - subtrahend.data.astype('float64', casting='safe')
                  for minuend, subtrahend in zip(minuend_hdul, subtrahend_hdul)
                  if minuend.data is not None and subtrahend.data is not None]

    return difference


# filename1 = '/home/lee/Documents/bias_frames/c4d_170331_202028_zri.fits.fz'
# filename2 = '/home/lee/Documents/bias_frames/c4d_170331_191113_zri.fits.fz'

filename1 = '/home/lee/Documents/k4m_160531_050920_ori.fits.fz'
filename2 = '/home/lee/Documents/k4m_161228_132947_dri.fits.fz'

# Open a new ds9 instance, or if already open, access it
display = pyds9.DS9(target='display', start='-title display')

with fits.open(filename1) as minuend_hdul, fits.open(filename2) as subtrahend_hdul:
    # minuend_hdul.info()
    # subtrahend_hdul.info()

    data_list = frame_subtract(minuend_hdul, subtrahend_hdul)

    # for data in data_list:
    #     display.set('frame new')
    #     display.set_np2arr(data)


