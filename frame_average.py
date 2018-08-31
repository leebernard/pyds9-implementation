"""
A function for pixel by pixel averaging of an arbitrary number of frames
"""

import numpy as np
import os
import pyds9
import re
import datetime
import gc

from astropy.io import fits
from contextlib import ExitStack
from astropy.stats import sigma_clip
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats

from get_filenames import get_filenames
# import timing


def _frame_mean(hdul_list):
    # unpack the data
    image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    # generate a list that contains how many data points were used to produce the average in each pixel
    # this can be pregenerated in this case, because no data points will be rejected
    pixel_counter = [[np.ones(image.shape) * len(image_stack) for image in image_list] for image_list in image_stack]

    # transpose the lists of data, then take the mean
    frame_mean = [np.mean(np.stack([hdu_data for hdu_data in data_tuple ]), axis=0) for data_tuple in zip(*image_stack)]

    return frame_mean, pixel_counter


def _frame_median(hdul_list):
    # unpack the data
    image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    # transpose the lists, then stack and take the median
    return [np.median(np.stack([hdu_data for hdu_data in data_tuple]), axis=0) for data_tuple in zip(*image_stack)]


def _sigma_clipped_frame_average(hdul_list, **kwargs):
    # generate an array of zeros the size of the array
    image_values = []
    with hdul_list[0] as hdul:
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

    # unpack the data, filtering out None values
    image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    # generate an object for sigmaclipping. Hopefully, this is faster than function calling each time
    sigclip = SigmaClip(**kwargs)
    # iterate through the list of files, summing up all values of corresponding pixels
    for image_list in image_stack:
        for n, image_data in enumerate(image_list):
            # extract the sigma clipped data

            clipped_data = sigclip(image_data)
            # add the sigma clipped data to the image array average, with clipped values set to zero
            # uses tracker to align entries
            image_values[n] = image_values[n] + clipped_data.filled(0)

            # keep tract of which pixels had a value added to them
            pixel_counter[n] = pixel_counter[n] + ~clipped_data.mask

            # # test section: comment out or remove in final version
            # image_values[n - backstep] = image_values[n - backstep] + hdu.data
            # pixel_counter[n - backstep] = pixel_counter[n - backstep] + np.ones(image_values[n - backstep].shape)
            # display.set_np2arr(hdu.data)

    # make a copy, to preserve how many pixels got counted zero times.
    # The pixels with zero counts are pixels that got sigma clipped out entirely
    average_counter = pixel_counter.copy()
    # set zero values to 1, to prevent dividing zero by zero
    for counter in average_counter:
        counter[counter == 0] = 1

    # compute average
    average_image = [image / counter for image, counter in zip(image_values, average_counter)]

    return average_image, pixel_counter


def frame_average(filename_list, sigma_clip=False, **kwargs):
    """
    This function returns an average of images stored in fits files, by
    accessing them via a list of filenames.

    This function averages frame data that has been stored on disk in fits
    files. It does so by taking a list of file names. If the files are not in
    the current directory, the file path should be included in the file name.
    This function presumes that each fits file contains a Header Data Unit
    with extensions: the averages are computed on a per extension basis, and
    returned as a list of data arrays. How many data points are used to
    calculate the average is kept track of on a pixel basis, and returned as a
    list of integer arrays that match the size of the list of image data arrays.

    There is an option for sigma clipped: however, sigma clipping significantly
    increases the runtime. If a particular pixel is clipped out in all images,
    the average of that pixel will be returned as zero. It's corresponding
    location in the data points tracker will also be zero. If you wish to
    eliminate defects such as cosmic rays, but do not want to wait for the
    sigma clipping to finish, consider using frame_median.

    Parameters
    ----------
    filename_list: list of strings
        Strings containing the file names of the frames to be averaged.
    sigma_clip: bool, optional
        If True, the data will be sigma clipped
    kwargs: dict, optional
        keyword arguments to be passed to sigma clipping

    Returns
    -------
    frame_mean: list
        A list of ndarrays containing image data aftera averaging

    """
    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(fits_name)) for fits_name in filename_list]
        if sigma_clip:
            return _sigma_clipped_frame_average(hdul_list, **kwargs)
        else:
            return _frame_mean(hdul_list)

def frame_median(filename_list):
    """
    This function returns the median of several frames.

    The the frames are access from disk, by passing the file names. The median
    is intended to be a fast way of getting a bias frame, while removing any
    outliers like cosmic rays.

    Parameters
    ----------
    filename_list: list of strings
        Strings containing the file names of the frames to be averaged.

    Returns
    -------
        list of ndarrays: frame data on a per amplifier basis, after taking the median
    """

    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(fits_name)) for fits_name in fits_list]
        return _frame_median(hdul_list)
"""
test code:
list1 = [1, 2, 3, 4, 5]
list2 = ['a', 'b', 'c', 'd', 'e']
list3 = [6, 7, 8, 9, 10]
overlist = [list1, list2, list3]

newlist = [ [overlist[i] for i in range(len(list1)]
print(newlist)
newlist = [[entry for entry in test_tuple] for test_tuple in zip(*overlist)]
print(newlist)
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

# a better way of opening multiple files
with ExitStack() as fits_stack:
    hdul_list = [fits_stack.enter_context(fits.open(fits_name, memmap=False)) for fits_name in fits_list]

    # image_average = _sigma_clipped_frame_average(hdul_list)
    # image_median = _frame_median(hdul_list)
    image_average = _frame_mean(hdul_list)

"""
# frame_average, frame_counts = sigma_clipped_frame_average(fits_list, sigma=5)

# try getting the median instead, by stacking the arrays
frame_average = _frame_median(fits_list)
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
        hdu.header.add_comment(str(datetime.date.today()))
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
"""
