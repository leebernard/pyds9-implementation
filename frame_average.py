"""
A function for pixel by pixel averaging of an arbitrary number of frames
"""

import numpy as np
import os
import pyds9

from astropy.io import fits
from astropy.stats import sigma_clip
biasframe_path = '/home/lee/Documents/bias_frames'

filename_list = os.listdir(biasframe_path)

# ensure all filenames have the proper extension
fits_list = [biasframe_path + '/' + filename for filename in filename_list if filename[-8:] in ['.fits.fz']]

# generate an array of zeros the size of the array
average_image = np.zeros(fits.getdata(fits_list[0]).shape)
# generate an array that keeps track of how many pixels have been added, on a pixel by pixel basis
pixel_counter = np.copy(average_image)

# iterate through the list of files, summing up all values of corresponding pixels
for fits_file in fits_list:
    with fits.open(fits_file) as hdul:
        # extract the sigma clipped data
        print(hdul[1].data.shape)
        clipped_data = sigma_clip(hdul[2].data, sigma=5)
        # add the sigma clipped data to the image array average, with clipped values set to zero
        average_image = average_image + clipped_data.filled(0)
        # keep tract of which pixels had a value added to them
        pixel_counter = pixel_counter + ~clipped_data.mask

# set zero values to 1, to prevent dividing zero by zero
pixel_counter[pixel_counter == 0] = 1

# compute average
return_value = average_image/pixel_counter

ds9 = pyds9.DS9()

print(np.histogram(pixel_counter, bins=np.arange(25)+.5))
ds9.set_np2arr(return_value)
