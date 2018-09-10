"""
A function for pixel by pixel averaging of an arbitrary number of frames
"""

import numpy as np
# import os
import pyds9
# import re
import datetime
import gc

from astropy.io import fits
from contextlib import ExitStack
from astropy.stats import sigma_clip
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats

from get_filenames import get_filenames
import timing


def _frame_mean(image_stack):
    # transpose the lists of data, then take the mean
    frame_mean = [np.mean(np.stack([hdu_data for hdu_data in data_tuple]), axis=0) for data_tuple in zip(*image_stack)]
    # transpose the lists of data, than take the std deviation and calculate error.
    # the error is calculated as std_dev/sqrt(n)
    pixel_error = [np.std(np.stack([hdu_data for hdu_data in data_tuple]) / np.sqrt(len(image_stack)), axis=0) for
                   data_tuple in zip(*image_stack)]

    return frame_mean, pixel_error


def _frame_median(image_stack):

    # transpose the lists, then stack and take the median
    frame_median = [np.median(np.stack([hdu_data for hdu_data in data_tuple]), axis=0) for data_tuple in zip(*image_stack)]

    # transpose the lists of data, than take the std deviation.
    # the error is calculated as 1.253*std_dev/sqrt(n)
    # Note that this is correct only for normal distributions, which is almost certainly not the case
    # Therefore these error values are suspect
    pixel_error = [1.253 * np.std(np.stack([hdu_data for hdu_data in data_tuple]) / np.sqrt(len(image_stack)), axis=0)
                   for data_tuple in zip(*image_stack)]

    return frame_median, pixel_error


def _sigma_clipped_frame_average(image_stack, **kwargs):
    # generate an array of zeros the size of the array
    image_values = [np.zeros_like(image_data, dtype=float) for image_data in image_stack[0]]

    variance_tracker = np.zeros(len(image_stack[0]))
    pixel_counter = [np.zeros_like(image_data, dtype=int) for image_data in image_stack[0]]

    # generate an object for sigmaclipping. Hopefully, this is faster than function calling each time
    sigclip = SigmaClip(**kwargs)
    # iterate through the list of files, summing up all values of corresponding pixels
    for image_list in image_stack:
        print('Variance:', variance_tracker)
        for n, image_data in enumerate(image_list):
            # extract the sigma clipped data

            clipped_data = sigclip(image_data)
            # add the sigma clipped data to the image array average, with clipped values set to zero
            # uses tracker to align entries
            image_values[n] = image_values[n] + clipped_data.filled(0)

            # keep track of the average variance
            variance_tracker[n] += np.ma.var(clipped_data)/len(image_stack)

            # keep tract of which pixels had a value added to them
            pixel_counter[n] = pixel_counter[n] + ~clipped_data.mask

    # make a copy, to preserve how many pixels got counted zero times.
    # The pixels with zero counts are pixels that got sigma clipped out entirely
    average_counter = pixel_counter.copy()
    # set zero values to 1, to prevent dividing zero by zero
    for counter in average_counter:
        counter[counter == 0] = 1

    # compute average
    average_image = [image/counter for image, counter in zip(image_values, average_counter)]

    # compute error
    pixel_error = [np.sqrt(image_variance/counter) for image_variance, counter in zip(variance_tracker, average_counter)]

    return average_image, pixel_error, pixel_counter


def _copy_hdul(hdul):
    # iterate through the HDUList, copying each header data unit
    hdu_list = [hdu.copy() for hdu in hdul]
    # generate and return a new HDUList
    return fits.HDUList(hdu_list)


def _write_average_data_to_file(data_list, writeto_filename, source_file_list=None, file_path='.', comment_string=None):
    """
    This is a specific method for frame averaging.

    This is a method specific to averaging frames, so this docstring is not
    published as part of the documentation. However, this method contains
    important techniques for creating new HDULs from data derived from other
    HDULs, so it is being given formal documentation.

    This function takes the Header Data Units from the entry of a list of file
    names, and replaces the data with the data list provided. It then modifies
    the Primary header to indicate that the data is an average of several
    frames, and adds the source filenames as a comment. It then writes the new
    HDUL to file, with the specified file name and directory.

    Parameters
    ----------
    data_list
    writeto_filename
    source_file_list
    file_path
    comment_string

    Returns
    -------

    """

    # add a None to the start of the data list, for the primary HDU
    data_list.insert(0, None)
    # copy the first HDUList
    with fits.open(file_path + '/' + source_file_list[0]) as hdul:
        hdulcopy = _copy_hdul(hdul)

    # make generator for modifying the fits file
    hdul_generator = (hdu for hdu in hdulcopy)

    for hdu, avg_data in zip(hdul_generator, data_list):
        hdu.data = avg_data
        hdu.header.set('OBSTYPE', 'average zero')
        # input('press Enter to continue')
        if type(hdu) is fits.hdu.image.PrimaryHDU:
            hdu.header.add_comment(str(datetime.date.today()))
            hdu.header.add_comment('Modified OBSTYPE')
            if comment_string:
                hdu.header.add_comment(comment_string)
            hdu.header.add_comment('Source file names:')
            for filename in filename_list:
                hdu.header.add_comment(filename)

    hdulcopy.writeto(file_path + '/' + writeto_filename)


def sigma_clipped_frame_average(filename_list, sigma=3.0, iters=5, **kwargs):
    """
    This function calculates a sigma clipped average of images stored in fits
    files, by accessing them via a list of filenames.

    This function averages frame data that has been stored on disk in fits
    files. It does so by taking a list of file names. If the files are not in
    the current directory, the file path should be included in the file name.
    This function presumes that each fits file contains a Header Data Unit
    with extensions: the averages are computed on a per extension basis, and
    returned as a list of data arrays.

    Pixels are clipped by comparing the pixel to the standard
    deviation of the image the pixel is in. How many data points are used to
    calculate the average is kept track of on a pixel basis, and is returned
    as a list of integer arrays that matches the list of image data arrays. If
    a particular pixel is clipped out in all images, the average of that pixel
    will be returned as zero. It's corresponding location in the data points
    tracker will also be zero.

    Error is calculated by taking the average variance of the clipped image
    data, and dividing the result by the pixel tracker. This produces an array
    of the same shape as the original image, with entries containing variances
    on a pixel by pixel basis that have been reduced by the number of pixels
    used to calculate the average for that pixel. The elementwise square root
    of this array is returned as the error.

    Sigma clipping significantly increases the runtime. If you wish to
    eliminate defects such as cosmic rays, but do not want to wait for the
    sigma clipping to finish, consider using frame_median.

     Parameters
    ----------
    filename_list: list of strings
        Strings containing the file names of the frames to be averaged.
    sigma: float, optional
        The number of standard deviations to use for both the lower and upper
        clipping limit.
    iters: int or None, optional
        The number of iterations to perform sigma clipping, or None to clip
        until convergence is achieved (i.e., continue until the last iteration
        clips nothing).
    kwargs: dict, optional
        keyword options for sigma clipping. See Astropy SigmaClip.

    Returns
    -------
    frame_mean: list
        A list of ndarrays containing image data after averaging
    pixel_error: list
        A list of ndarrays containing the error on the average of each pixel.
    pixel_tracker: list
        A list of ndarrays containing how many data point were used to produce
        an average, on a pixel by pixel basis.
    """
    # unpack the data, ignoring None values
    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    return _sigma_clipped_frame_average(image_stack, sigma, iters, **kwargs)


def frame_average(filename_list):
    """
    This function returns an average and error estimate of images stored in fits files, by
    accessing them via a list of filenames.

    This function averages frame data that has been stored on disk in fits
    files. It does so by taking a list of file names. If the files are not in
    the current directory, the file path should be included in the file name.
    This function presumes that each fits file contains a Header Data Unit
    with extensions. Data is taken from the extensions, stacked, and the
    averages computed. This produces a list of data arrays corresponding to
    the HDU data arrays.

    The error is calculated by taking the standard deviation of the data used
    to calculate the mean, on a pixel by pixel basis, and then dividing by the
    square root of the number of frames. This produces a list of arrays the
    same shape as the frame data.

    Parameters
    ----------
    filename_list: list of strings
        Strings containing the file names of the frames to be averaged.

    Returns
    -------
    frame_mean: list
        A list of ndarrays containing image data after averaging
    pixel_error: list
        A list of ndarrays correspoding to frame_mean, that contains the
        per-pixel error on the average.
    """
    # unpack the data, ignoring None values
    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    return_value = _frame_mean(image_stack)

    # clean up the mess
    print('Amount of garbage:')
    print(gc.collect())
    return return_value


def frame_median(filename_list, writeto_filename=None):
    """
    This function returns the median of several frames.

    The frames are accessed from disk, by passing the file names. If the files
    are not in the current directory, the file path should be included in the
    file name. This function presumes that each fits file contains a Header
    Data Unit with extensions. The data is stacked, and the median calculated
    on a per pixel basis. The median is intended to be a fast way of getting
    a bias frame, with some robustness against any outliers like cosmic rays.

    The error is calculated by finding the standard error on the mean, and
    multiplying by 1.253. Note that this only gives accurate error estimates
    for data that is normally distributed. This should be kept in mind when
    reviewing errors on data sets that contain outliers.

    Parameters
    ----------
    filename_list: list of strings
        Strings containing the file names of the frames to be averaged.

    Returns
    -------
        list of ndarrays: frame data on a per amplifier basis, after taking the median
    """

    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    return _frame_median(image_stack)


'''
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
'''

biasframe_path = '/home/lee/Documents/bias_frames'
extension = '.fits.fz'
pattern = '(?=.*k4m)'  # look-ahead regex pattern that checks for 'k4m'

# retrieve the filenames
filename_list = get_filenames(biasframe_path, extension=extension, pattern=pattern, include_path=False)

fits_list = [biasframe_path + '/' + filename for filename in filename_list]

image_average, pixel_error = frame_median(fits_list)
# image_average, pixel_error = frame_average(fits_list)
# image_average, pixel_error, pixel_tracker = frame_average(fits_list, sigma_clip=True, sigma=5)

# print('Amount of garbage:')
# print(gc.collect())

ds9 = pyds9.DS9(target='ds9')
for image in image_average:
    ds9.set('frame new')
    try:
        ds9.set_np2arr(image)
    except ValueError:
        print('invalid data')
        print('data:', image)

#     # image_average = _sigma_clipped_frame_average(hdul_list)
#     image_median = _frame_median(hdul_list)
#     # image_average = _frame_mean(hdul_list)




# a better way of opening multiple files
with ExitStack() as fits_stack:
    hdul_list = [fits_stack.enter_context(fits.open(fits_name)) for fits_name in fits_list]
    image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

# frame_average, frame_counts = sigma_clipped_frame_average(fits_list, sigma=5)

# # try getting the median instead, by stacking the arrays
# frame_average = _frame_median(fits_list)
# # insert a None, for the Primary HDU
# frame_average.insert(0, None)

# for displaying the raw data
# display = pyds9.DS9(target='display', start='-title display')


_write_average_data_to_file(average_data, 'testaverage.fits.fz', source_file_list=filename_list,
                            file_path=biasframe_path, comment_string='Changed data to the median of 15 zero frames, of which this is the first')

ds9 = pyds9.DS9(target='ds9')

for image in frame_average:
    ds9.set('frame new')
    ds9.set_np2arr(image)

