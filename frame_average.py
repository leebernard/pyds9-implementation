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
from astropy.stats import sigma_clipped_stats

from get_filenames import get_filenames
# import timing


def frame_mean(filename_list):
    """
    Returns the mean of a stack of images, loaded from disk.


    Parameters
    ----------
    filename_list

    Returns
    -------

    """
    # make a list of empty lists to hold the data
    with fits.open(filename_list[0]) as hdul:
        image_stacks = [[] for hdu in hdul if hdu.data is not None]

    # unpack the data in a fashion that can be stacked
    for fits_file in filename_list:
        with fits.open(fits_file) as hdul:
            backstep = 0
            for n, hdu in enumerate(hdul):
                if hdu.data is None:
                    backstep += 1
                else:
                    image_stacks[n - backstep].append(hdu.data)

    # image_stacks = [np.stack(stack) for stack in image_stacks]

    return [np.mean(np.stack(stack), axis=0) for stack in image_stacks]


def sigma_clipped_frame_average(filename_list, **kwargs):
    # generate an array of zeros the size of the array
    image_values = []
    with fits.open(filename_list[0]) as hdul:
        hdul.info()
        for hdu in hdul:
            if hdu.data is None:
                print('No Data found, skipped')
            else:
                print(hdu.data.shape)
                image_values.append(np.zeros(hdu.data.shape))
        # generate an array that keeps track of how many pixels have been added, on a pixel by pixel basis
    pixel_counter = image_values.copy()

    # for displaying the raw data
    # display = pyds9.DS9(target='display', start='-title display')

    # generate an object for sigmaclipping. Hopefully, this is faster than function calling each time
    sigclip = SigmaClip(**kwargs)
    # iterate through the list of files, summing up all values of corresponding pixels
    for fits_file in filename_list:
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

    # make a copy, to preserve how many pixels got counted zero times.
    # The pixels with zero counts are pixels that got sigma clipped out entirely
    average_counter = pixel_counter.copy()
    # set zero values to 1, to prevent dividing zero by zero
    for counter in average_counter:
        counter[counter == 0] = 1

    # compute average
    average_image = [image / counter for image, counter in zip(image_values, average_counter)]

    return average_image, pixel_counter


def frame_median(filename_list):
    # make a list of empty lists to hold the data
    with fits.open(filename_list[0]) as hdul:
        image_stacks = [[] for hdu in hdul if hdu.data is not None]

    # unpack the data in a fashion that can be stacked
    for fits_file in fits_list:
        with fits.open(fits_file) as hdul:
            backstep = 0
            for n, hdu in enumerate(hdul):
                if hdu.data is None:
                    backstep += 1
                else:
                    image_stacks[n - backstep].append(hdu.data)

    # image_stacks = [np.stack(stack) for stack in image_stacks]

    frame_median = [np.median(np.stack(stack), axis=0) for stack in image_stacks]

    return frame_median


"""
test code:
list1 = [1, 2, 3, 4, 5]
list2 = ['a', 'b', 'c', 'd', 'e']
list3 = [6, 7, 8, 9, 10]
overlist = [list1, list2, list3]
len(overlist)
3
len(overlist[0])
5
overlist[0]
[1, 2, 3, 4, 5]
overlist[:]
[[1, 2, 3, 4, 5], ['a', 'b', 'c', 'd', 'e'], [6, 7, 8, 9, 10]]
overlist[:][0]
[1, 2, 3, 4, 5]
overlist[0][:]
[1, 2, 3, 4, 5]
"""

biasframe_path = '/home/lee/Documents/bias_frames'
extension = '.fits.fz'
pattern = '(?=.*k4m)'

# retrieve the filenames
filename_list = get_filenames(biasframe_path, extension=extension, pattern=pattern, include_path=False)

fits_list = [biasframe_path + '/' + filename for filename in filename_list]

# frame_average, frame_counts = sigma_clipped_frame_average(fits_list, sigma=5)

# try getting the median instead, by stacking the arrays
frame_average = frame_median(fits_list)
# insert a None, for the Primary HDU
frame_average.insert(0, None)

# for displaying the raw data
# display = pyds9.DS9(target='display', start='-title display')


with fits.open(fits_list[0]) as hdul:
    # hdu_list = [hdu.copy() for hdu in hdul]
    hdulcopy = fits.HDUList([hdu.copy() for hdu in hdul])

# make generator for modyfying the fits file
hdul_generator = (hdu for hdu in hdulcopy)



for hdu, avg_data in zip(hdul_generator, frame_average):
    hdu.data = avg_data
    hdu.header.set('OBSTYPE', 'average zero')
    # input('press Enter to continue')
    if type(hdu) is fits.hdu.image.PrimaryHDU:
        hdu.header.add_comment('Modified OBSTYPE')
        hdu.header.add_comment('Changed data to the median of 15 zero frames, of which this is the first')
        hdu.header.add_comment('Source file names:')
        for filename in filename_list:
            hdu.header.add_comment(filename)

hdulcopy.writeto(biasframe_path + '/testaverage' + extension)

ds9 = pyds9.DS9(target='ds9')

for image in frame_average:
    ds9.set('frame new')
    ds9.set_np2arr(image)

