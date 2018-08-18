"""
A function for pixel by pixel averaging of an arbitrary number of frames
"""

import numpy as np
import os

from astropy.io import fits

biasframe_path = '/home/lee/Documents/bias_frames'

filename_list = os.listdir(biasframe_path)

# ensure all filenames have the proper extension
fits_list = [biasframe_path + '/' + filename for filename in filename_list if filename[-8:] in ['.fits.fz']]

# acquire the size of the image, number of tiles, etc
hdul = fits.open(fits_list[0])
