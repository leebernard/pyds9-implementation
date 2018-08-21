"""
A function for pixel by pixel averaging of an arbitrary number of frames
"""

import numpy as np
import os
import pyds9
import re

from astropy.io import fits
from astropy.stats import sigma_clip


def get_filenames(path, extension='', pattern=''):
    # retrieve all filenames from the directory
    filename_list = os.listdir(biasframe_path)

    # convert extension and pattern to raw strings
    r_extension = "%r"%extension
    # ensure all filenames have the proper extension

    fits_list = [biasframe_path + '/' + filename for filename in filename_list if
                 re.search(extension+r'$', filename) and re.search(pattern, filename)]

    return fits_list


biasframe_path = '/home/lee/Documents/bias_frames'
extension = r'\.fits\.fz'
pattern = '(?=.*k4m)'

# retrieve the filenames
fits_list = get_filenames(biasframe_path, extension=extension, pattern=pattern)

# generate an array of zeros the size of the array
average_image = []
with fits.open(fits_list[0]) as hdul:
    hdul.info()
    for hdu in hdul:
        if type(hdu.data) is np.ndarray:
            print(hdu.data.shape)
            average_image.append(np.zeros(hdu.data.shape))
    # generate an array that keeps track of how many pixels have been added, on a pixel by pixel basis
pixel_counter = average_image.copy()


# for displaying the raw data
display = pyds9.DS9(target='display', start='-title display')
# iterate through the list of files, summing up all values of corresponding pixels
for fits_file in fits_list:
    with fits.open(fits_file) as hdul:

        for hdu in hdul:

            # extract the sigma clipped data
            clipped_data = sigma_clip(hdu.data, sigma=5)
            # add the sigma clipped data to the image array average, with clipped values set to zero
            average_image = average_image + clipped_data.filled(0)

            # keep tract of which pixels had a value added to them
            pixel_counter = pixel_counter + ~clipped_data.mask

            # test section: comment out or remove in final version
            # average_image = average_image + hdu.data
            # pixel_counter = pixel_counter + np.ones(average_image.shape)
            # display.set_np2arr(hdu.data)
            # print(hdu.data.shape)

print(np.histogram(pixel_counter, bins=np.arange(25)))
# set zero values to 1, to prevent dividing zero by zero
pixel_counter[pixel_counter == 0] = 1

# compute average
return_value = average_image/pixel_counter

ds9 = pyds9.DS9(target='ds9')


ds9.set('frame new')
ds9.set_np2arr(return_value)


