"""
A function for pixel by pixel averaging of an arbitrary number of frames
"""

import numpy as np
import os
import pyds9
import re

from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.stats import SigmaClip

import timing

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
image_values = []
with fits.open(fits_list[0]) as hdul:
    hdul.info()
    for hdu in hdul:
        if hdu.data is None:
            print('Data is', hdu.data, ', skipped')
        else:
            print(hdu.data.shape)
            image_values.append(np.zeros(hdu.data.shape))
    # generate an array that keeps track of how many pixels have been added, on a pixel by pixel basis
pixel_counter = image_values.copy()


# for displaying the raw data
# display = pyds9.DS9(target='display', start='-title display')

# generate an object for sigmaclipping. Hopefully, this is faster than function calling each time
sigclip = SigmaClip(sigma=5, iters=5)
# iterate through the list of files, summing up all values of corresponding pixels
for fits_file in fits_list:
    with fits.open(fits_file) as hdul:
        # a way of keeping track of how many hdu's have been skipped due to having no data
        backstep = 0
        for n, hdu in enumerate(hdul):
            print(backstep)
            if hdu.data is None:
                # if there is no data entry for this hdu, skip it while keeping track
                backstep += 1
            else:
                # extract the sigma clipped data

                clipped_data = sigclip(hdu.data)
                # add the sigma clipped data to the image array average, with clipped values set to zero
                # uses tracker to align entries
                image_values[n - backstep] = image_values[n - backstep] + clipped_data.filled(0)

                # keep tract of which pixels had a value added to them
                pixel_counter[n - backstep] = pixel_counter[n - backstep] + ~clipped_data.mask

                # # test section: comment out or remove in final version
                # image_values[n - backstep] = image_values[n - backstep] + hdu.data
                # pixel_counter[n - backstep] = pixel_counter[n - backstep] + np.ones(image_values[n - backstep].shape)
                # display.set_np2arr(hdu.data)
        print('hdu:', hdul.filename())

print(np.histogram(pixel_counter, bins=np.arange(25)))
# set zero values to 1, to prevent dividing zero by zero
for counter in pixel_counter:
    counter[counter == 0] = 1

# compute average
average_image = [image/counter for image, counter in zip(image_values, pixel_counter)]


ds9 = pyds9.DS9(target='ds9')

for image in average_image:
    ds9.set('frame new')
    ds9.set_np2arr(image)


