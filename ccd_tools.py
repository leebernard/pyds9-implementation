"""
''ccd_tools'' is a wrapper and extension for ''pyds9'', with basic stats function
for convenience. This is very much still in development, so use with
caution.

Classes
-------
The 'Region' class is for holding segments of image data in a convenient
wrapper.

Functions
---------
get_ds9_regions:
    Retrieves a selected region from DS9, and returns it as a Region instance
make_source_mask:
    Generates and returns a source mask from image data.
image_stats:
    Returns stats on image data.
bias_from_ds9:
    Returns the data in the bias section from an image loaded in DS9.
sky_subtract:
    Returns the sky background subtracted image data
bias_subtract:
    Returns the bias subtracted data
display_data:
    A wrapper for displaying image data using MatPlotLib
get_filenames:
    A function that retrieves file names from a particular directory.
frame_subtract:
    Subtracts image data frame by frame, either from file or from DS9
frame_average:
    Takes the average of a set of image data, from file.
frame_median:
    Takes the per-pixel median of a set of image data, from file. Intended
    to be a more robust version of frame_average.
sigma_clipped_frame_average:
    Takes the average of a set of image data from file, with outliers
    clipped.
"""

__version__ = '0.2'
__author__ = 'Lee Bernard'
# import needed packages

import matplotlib.pyplot as plt
import warnings
import numpy as np
import re
import datetime
import gc
import os

from contextlib import ExitStack
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
# from astropy.stats import sigma_clip
from astropy.stats import SigmaClip
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

import pyds9


class Region:
    """
    This class is for convenient packaging of the region data.

    For example, (self.xmin, self.ymin) gives the array location of lower left corner of region in the image data.

    Attributes
    ----------
    region_def: string
        String definition of the region. This can be the region as defined in DS9, or any string definition the user
        decides.
    source_file: string
        String containing the name and path of the file that the region originated from.
    x_coord: float
        x coordinate of the region center
    y_coord: float
        y coordinate of the region center
    width: float
        width of the region
    height: float
        height of the region
    xmin: int
        x coordinate of the lower left corner
    xmax: int
        x coordinate of the upper right corner
    ymin: int
        y coordinate of the lower left corner
    ymax: int
        y coordinate of the upper right corner
    data: ndarray
        Image data array of the defined region.
    bias_sub: ndarray
        Image data array of the defined region, after bias subtraction
    sky_sub: ndarray
        Image data array of the defined region, after background subtraction.
        This data may also be bias subtracted.
    bias_stats: tuple
        Three numbers representing, in order: the bias region mean, median, and standard deviation

    Methods
    -------
    sky_subtract: adds sky background subtracted data to the sky_sub attribute.
    stats: displays statistics on the data attributes
    """

    def __init__(self):
        # set data members so that if any are not being used they return none

        # important meta data
        self.region_def = None
        self.source_file = None

        self.x_coord = None
        self.y_coord = None
        self.width = None
        self.height = None

        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None

        self.data = None
        self.bias_sub = None
        self.sky_sub = None

        self.bias_stats = None
        self.sky_stats = None

    def stats(self, mask=None, sigma_clip=False, mask_sources=False,  **kwargs):
        """
        Displays statistics of the region data, and bias subtracted and sky
        subtracted data if available.

        Has options for masking and sigma clipping.

        Parameters
        ----------
        sigma_clip: bool, optional
            If True, returns sigma clipped statistics.
        mask: np.ndarray (bool), optional
            Boolean array of the same shape as self.data. Entries that are set to
            True are ignored in statistical calculations.
        mask_sources:
            Flag for masking sources or not
        kwargs: dict
            Keyword arguments to be passed to image stats

        See Also
        --------
        image_stats: calculates statistics on image data
        """

        print('-----------------------')
        print('Region Data Statistics:')
        image_stats(self.data, mask=mask, sigma_clip=sigma_clip, mask_sources=mask_sources, **kwargs)
        # if there are stats on the bias and sky background, print those also
        if isinstance(self.bias_sub, np.ndarray):

            print('---------------------')
            print('Region Bias Subtracted Statistics:')
            image_stats(self.bias_sub, mask=mask, sigma_clip=sigma_clip, mask_sources=mask_sources, **kwargs)
        if isinstance(self.sky_sub, np.ndarray):

            print('---------------------')
            print('Sky Subtracted Data Statistics:')
            image_stats(self.sky_sub, mask=mask, sigma_clip=sigma_clip, mask_sources=mask_sources, **kwargs)

        # return mean, median, std

    def sky_subtract(self, mask=None, **kwargs):
        """
        This method calculates the background subtracted data, and stores it in
        the sky_sub attribute.

        This is a wrapper for sky_subtract. If bias subtracted data is available,
        it uses that. Otherwise, it uses data in the data attribute. If neither
        are available, raises an exception.

        Parameters
        ----------
        mask: numpy.ndarray (bool), optional
            A boolean mask with the same shape as data, where a True value
            indicates the corresponding element of data is masked. Masked pixels
            are excluded when computing the statistics.
        kwargs: dict, optional
            Optional keywords for sky_subtract

        Raises
        ------
        ValueError
            If both the data attribute and the bias_sub attribute are set to none.

        See Also
        --------
        sky_subtract: returns background subtracted data
        """

        # first look for bias subtracted data
        if isinstance(self.bias_sub, np.ndarray):
            self.sky_sub, self.sky_stats = sky_subtract(self.bias_sub, mask=mask, **kwargs)
        # if bias subtracted data is not found, used raw data
        elif isinstance(self.data, np.ndarray):
            warnings.warn('Bias subtracted data not found. Using raw data.', category=UserWarning)
            self.sky_sub, self.sky_stats = sky_subtract(self.data, mask=mask, **kwargs)
        # if no data is available, throw an exception
        elif self.data is None:
            raise ValueError('Data attributes are unassigned')
        else:
            raise TypeError('Invalid value stored in data attributes.')


def _find_hdu_extension(extname, hdulist):
    # finds the first match in the hdulist for the extension name
    # This is needed for frame_subtracting using ds9, but should work for any format of hdu list
    for hdu in hdulist:
        try:
            if re.match(extname, hdu.header['extname'], flags=re.IGNORECASE):
                hdu_match = hdu
                # once the extension is found, break
                break
        except KeyError:
            # catch the exception when no extension name exists, and keep moving
            print('No extension name.')
    else:
        # corresponding extension not found, raise error
        raise RuntimeError('A matching extension was not found.')

    return hdu_match


def _frame_subtract(minuend_hdul, subtrahend_hdul):
    return [minuend.data.astype('float64', casting='safe') - subtrahend.data.astype('float64', casting='safe')
            for minuend, subtrahend in zip(minuend_hdul, subtrahend_hdul)
            if minuend.data is not None and subtrahend.data is not None]


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
    image_median_list = [np.median(np.stack([hdu_data for hdu_data in data_tuple]), axis=0) for data_tuple in zip(*image_stack)]

    # transpose the lists of data, than take the std deviation.
    # the error is calculated as 1.253*std_dev/sqrt(n)
    # Note that this is correct only for normal distributions, which is almost certainly not the case
    # Therefore these error values are suspect
    pixel_error_list = [1.253 * np.std(np.stack([hdu_data for hdu_data in data_tuple]) / np.sqrt(len(image_stack)), axis=0)
                   for data_tuple in zip(*image_stack)]

    return image_median_list, pixel_error_list


def _sigma_clipped_frame_average(image_stack, **kwargs):
    # generate an array of zeros the size of the array
    image_value_sum_list = [np.zeros_like(image_data, dtype=float) for image_data in image_stack[0]]

    variance_tracker_list = np.zeros(len(image_stack[0]))
    pixel_counter_list = [np.zeros_like(image_data, dtype=int) for image_data in image_stack[0]]

    # generate an object for sigmaclipping. Hopefully, this is faster than function calling each time
    sigclip = SigmaClip(**kwargs)
    # iterate through the list of files, summing up all values of corresponding pixels
    for image_list in image_stack:
        # # test code for checking calculations. Comment out or remove in final version
        # print('Variance:', variance_tracker)
        for image_data, image_value_sum, variance_tracker, pixel_counter \
                in zip(image_list, image_value_sum_list, variance_tracker_list, pixel_counter_list):
            # extract the sigma clipped data

            clipped_data = sigclip(image_data)
            # add the sigma clipped data to the image array average, with clipped values set to zero
            # uses tracker to align entries
            image_value_sum += clipped_data.filled(0)

            # keep track of the average variance
            variance_tracker += np.ma.var(clipped_data)/len(image_stack)

            # keep tract of which pixels had a value added to them
            pixel_counter += ~clipped_data.mask

    # make a copy, to preserve how many pixels got counted zero times.
    # The pixels with zero counts are pixels that got sigma clipped out entirely
    average_counter_list = pixel_counter_list.copy()
    # set zero values to 1, to prevent dividing zero by zero
    for counter in average_counter_list:
        counter[counter == 0] = 1

    # compute average
    average_image_list = [image/counter for image, counter in zip(image_value_sum_list, average_counter_list)]

    # compute error
    pixel_error_list = [np.sqrt(image_variance/counter) for image_variance, counter in zip(variance_tracker_list, average_counter_list)]

    return average_image_list, pixel_error_list, pixel_counter_list


def _copy_hdul(hdul):
    # iterate through the HDUList, copying each header data unit
    hdu_list = [hdu.copy() for hdu in hdul]
    # generate and return a new HDUList
    return fits.HDUList(hdu_list)


def _write_average_data_to_file(data_list, writeto_filename, source_filename_list=None, file_path='.', comment_string=None):
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
    source_filename_list
    file_path
    comment_string

    Returns
    -------

    """

    # add a None to the start of the data list, for the primary HDU
    data_list.insert(0, None)
    # copy the first HDUList
    with fits.open(file_path + '/' + source_filename_list[0]) as hdul:
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
            for filename in source_filename_list:
                hdu.header.add_comment(filename)

    hdulcopy.writeto(file_path + '/' + writeto_filename)


def get_ds9_region(ds9=None, bias_sec=None, get_data=True):
    """
    Gets the first single valid box region selected in ds9, and returns it as
    a Region object.

    This function requires that the region be selected and defined as a box.
    (see DS9) It stores the region data and the bias subtracted data in a
    region instance, as well as the region definition, source file, and
    dimensions of the region.

    Parameters
    ----------
    ds9: DS9 object, optional
        optional parameter to specify a DS9 target
    get_data: bool, optional
        If True, does not retrieve data from DS9. This to reduce the resource
        requirements if data is being handled separately.
    bias_sec: array-like, optional
        This must be a list, tuple, or array of four numbers that define the
        bias section of the image, in the form (xmin, xmax, ymin, ymax). The
        numbers can be ints, floats, or strings that represent ints. If a float
        is given, it will be truncated towards zero. If None, the bias section
        will be retrieved from the header.

    Returns
    -------
    region: Region object

    Raises
    ------
    IndexError
        If a region has not been selected in DS9, or the region selected is not
        a box.
    Examples
    --------
    >>> from ccd_tools import *
    >>> ds9 = pyds9.DS9()
    >>> testregion = get_ds9_region(ds9=ds9)
    # Region file format: DS9 version 4.1
    global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    image
    Region definition:  ['box(893.07248,1213.0923,72.575962,72.575909,359.99999)']
    Bias from  /home/lee/Documents/k4m_160531_050920_ori.fits.fz[im1]
    Bias Section is [2066:2112,18:2065]
    ----------------
    Bias statistics:
    Mean: 3357.89
    Median: 3358.00
    Std: 4.50
    Bias statistics after bias subtraction:
    Mean: 0.00
    Median: 0.11
    Std: 4.50
    >>> testregion.data
    array([[5504, 5487, 5503, ..., 5518, 5527, 5563],
           [5488, 5495, 5527, ..., 5515, 5542, 5522],
           [5561, 5540, 5586, ..., 5514, 5558, 5552],
           ...,
           [5588, 5552, 5615, ..., 5548, 5579, 5545],
           [5534, 5530, 5553, ..., 5545, 5556, 5590],
           [5587, 5560, 5579, ..., 5526, 5560, 5519]], dtype=int32)
    >>> testregion.bias_sub
    array([[2146.10897334, 2129.10897334, 2145.10897334, ..., 2160.10897334,
            2169.10897334, 2205.10897334],
           [2130.10897334, 2137.10897334, 2169.10897334, ..., 2157.10897334,
            2184.10897334, 2164.10897334],
           [2203.10897334, 2182.10897334, 2228.10897334, ..., 2156.10897334,
            2200.10897334, 2194.10897334],
           ...,
           [2230.10897334, 2194.10897334, 2257.10897334, ..., 2190.10897334,
            2221.10897334, 2187.10897334],
           [2176.10897334, 2172.10897334, 2195.10897334, ..., 2187.10897334,
            2198.10897334, 2232.10897334],
           [2229.10897334, 2202.10897334, 2221.10897334, ..., 2168.10897334,
            2202.10897334, 2161.10897334]])


    """

    # if a ds9 target is not specified, make one
    if ds9 is None:
        try:
            ds9 = pyds9.DS9()
        except ValueError:
            print('Specify a target DS9() instance to retrieve the region from. '
                  'e.g:\n  d = pyds9.DS9(\'7f000101:43123\')\n  r = get_ds9_region(ds9=d)')
            raise

    # set the region format to ds9 default, and coordinate system to image. This ensures the format is standardized.
    # image format is required to properly index the data array.
    ds9.set('regions format ds9')
    ds9.set('regions system image')

    # get selected region info
    raw_string = ds9.get('regions selected')
    # print(raw_string)

    # transform string into list that is organized by lines
    pattern = re.compile('.+')  # regex pattern that finds everything not a newline
    str_list = pattern.findall(raw_string)

    try:
        while not re.match('box', str_list[0]):
            # if re.match('# tile', str_list[0]):
            #     print('Tile mode detected')
            #     tiled = True

            print(str_list.pop(0))
    # if the loop runs through the whole string without finding a region, print a message
    except IndexError:
        message = 'No valid region found. Please select a valid box region.'
        print(message)
        raise

    print('Region definition: ', str_list)
    # parse the meta data string
    # pattern finds all sequences of digits that may or may not contain a period
    pattern = re.compile('\d+\.?\d*')
    region_def = pattern.findall(str_list[0])  # should probably add an exception test here

    # make a region object to hold all the data
    region = Region()

    # region definition: origin is lower left, given as x and y coord, with a width and a height
    x_coord = float(region_def[0])
    y_coord = float(region_def[1])
    width = float(region_def[2])
    height = float(region_def[3])

    # region slicing data
    xmin = int(x_coord - width / 2)
    xmax = int(x_coord + width / 2)

    ymin = int(y_coord - height / 2)
    ymax = int(y_coord + height / 2)

    # if get_data is true, retrieve the data and bias subtracted data
    if get_data:

        with ds9.get_pyfits() as hdulist:
            hdu = hdulist[0]
            frame_data = hdu.data
            # determine what the valid data region is
            pattern = re.compile('\d+')  # pattern for all decimal digits
            try:
                datasec = pattern.findall(hdu.header['DATASEC'])
            except KeyError as err:
                message = str(err) + ' Unable to validate region.\nRegion may exceed image bounds'
                warnings.warn(message, category=UserWarning)
            else:
                # print('Data section: ', datasec)
                # if the selected region is outside the valid data range, throw a warning
                if int(datasec[0]) > xmin or int(datasec[1]) < xmax or int(datasec[2]) > ymin or int(datasec[3]) < ymax:
                    message = 'Region definition exceeds valid data range. \n ' \
                              'This could be due to entering the bias region, or crossing frames.'
                    warnings.warn(message, category=UserWarning)

        # store the raw data
        region.data = frame_data[ymin:ymax, xmin:xmax]

        # retrieve bias section, and calculate stats
        try:
            bias_data = bias_from_ds9(ds9, bias_sec)
        except KeyError as err:
            message = str(err) + ' Unable to perform bias subtraction.'
            warnings.warn(message, category=UserWarning)
        else:
            # bias dependent operations
            bias_mean, bias_median, bias_std = sigma_clipped_stats(bias_data, sigma=3.0, iters=5)
            print('----------------')
            print('Bias statistics:')
            print(f'Mean: {bias_mean:.2f}')
            print(f'Median: {bias_median:.2f}')
            print(f'Std: {bias_std:.2f}')
            # explicitly use the bias mean as the value to subtract
            bias = bias_mean
            # store the bias stats
            region.bias_stats = bias_mean, bias_median, bias_std

            # display the bias stats after subtraction
            test_mean, test_median, test_std = sigma_clipped_stats(bias_data - bias, sigma=3.0, iters=5)
            print('Bias statistics after bias subtraction: ')
            print(f'Mean: {test_mean:.2f}')
            print(f'Median: {test_median:.2f}')
            print(f'Std: {test_std:.2f}')

            # perform the bias subtraction
            bias_sub_frame = frame_data - bias
            region.bias_sub = bias_sub_frame[ymin:ymax, xmin:xmax]

    # package all the meta data
    region.x_coord = x_coord
    region.y_coord = y_coord
    region.width = width
    region.height = height

    region.xmin = xmin
    region.xmax = xmax
    region.ymin = ymin
    region.ymax = ymax

    region.source_file = ds9.get('file')
    region.region_def = str_list[0]

    return region


def make_source_mask(indata, snr=2, npixels=5, display_mask=False, **kwargs):
    """
    Wrapper for the make_source_mask from photutils. Makes a mask corresponding
    to any flux sources in the data.

    make_source_mask is imported as part of the function call, so if the
    photutils package is not installed, or broken, it does not break the rest
    of this package.

    Parameters
    ----------
    indata: array-like
        2D image data array containing a flux source to be masked.
    snr: float, optional
        Signal to noise ratio threshold above background to consider a pixel as
        being part of a source.
    npixels: int, optional
        Minimum number of continuous pixels that are above the threshold that
        an object must have to be considered a source.
    display_mask: bool, optional
        If True, a figure displaying the generated mask is produced
    kwargs: dict
        Keyword arguments for the make_source_mask function. See documentation
        for photutils.

    Returns
    -------
    bool_mask: ndarray, bool
        A boolean mask that corresponds to the data input. Pixels that are part
        of a source are set to True.
    """
    from photutils import make_source_mask

    mask = make_source_mask(indata, snr, npixels, **kwargs)

    if display_mask is True:
        plt.figure()
        plt.imshow(mask, origin='lower', cmap='viridis')
    return mask


def image_stats(imdata, mask=None, sigma_clip=False, mask_sources=False, **kwargs):
    """
    Calculates statistics on the background of the image data given.

    Parameters
    ----------
    imdata: array-like
        Data array of an image
    mask: numpy.ndarray (bool), optional
        A boolean mask with the same shape as data, where a True value
        indicates the corresponding element of data is masked. Masked pixels
        are excluded when computing the statistics.
    mask_sources: bool, optional
        If True, a source mask will be automatically generated using photutils.
        This source mask will be combined with any mask passed through the mask
        parameter. To set options for mask generation, generate a mask
        separately, then pass it through the mask parameter.
    sigma_clip: bool, optional
        If True, sigma clipped statistics will be returned using the function
        sigma_clipped_stats from astropy.stats
    kwargs: dict, optional
        Keyword arguments to be passed to sigma_clipped_stats

    Returns
    -------
    mean, median, stddev: float
        the mean, median, and standard deviation of the background
    """
    # if no mask is specified, mask any sources
    if mask_sources:
        obj_mask = make_source_mask(imdata, display_mask=True)
        # if both a mask is specified and object masking is called for, combine the masks
        if mask_sources and mask:
            mask = mask + obj_mask
        # otherwise, mask is just the object mask
        else:
            mask = obj_mask

    # apply the mask. If the mask is None, this effectively does nothing
    masked_data = np.ma.array(imdata, mask=mask)
    if sigma_clip:
        mean, median, std = sigma_clipped_stats(masked_data, **kwargs)
    else:
        mean, median, std = np.ma.mean(masked_data), np.ma.median(masked_data), np.ma.std(masked_data)

    print(f'Min pixel value: {np.min(masked_data):.2f}')
    print(f'Max pixel value: {np.max(masked_data):.2f}')
    print(f'Mean: {mean:.2f}')
    print(f'Median: {median:.2f}')
    print(f'Std: {std:.2f}')
    return mean, median, std


def bias_from_ds9(ds9=None, bias_sec=None):
    """
    This function retrieves the bias section from a fits file loaded in DS9.

    Parameters
    ----------
    ds9: DS9() instance, optional
        Target DS9 instance to retrieve the data from. If not specified, this
        will be acquired automatcially.
    bias_sec: array-like, optional
        This must be a list, tuple, or array of four numbers that define the
        bias section of the image, in the form (xmin, xmax, ymin, ymax). The
        numbers can be ints, floats, or strings that represent ints. If a float
        is given, it will be truncated towards zero. If None, the bias section
        will be retrieved from the header.

    Returns
    -------
    bias_data: ndarray
        Numpy array containing the pixel values of the bias section.
    """

    # if a ds9 target is not specified, make one
    if ds9 is None:
        try:
            ds9 = pyds9.DS9()
        except ValueError:
            print('Specify a target DS9() instance to retrieve the bias from. '
                  'e.g:\n  d = pyds9.DS9(\'7f000101:43123\')\n  r = get_ds9_region(ds9=d)')
            raise
    hdulist = ds9.get_pyfits()

    # hdulist.info()
    # extract header data unit from list
    hdu = hdulist[0]

    if bias_sec is None:
        # extract the bias definition
        bias_str = hdu.header['BIASSEC']

        print('Bias from ', ds9.get('file'))
        print('Bias Section is ' + bias_str)
        # print(type(bias_str))
        # slice the string, for converting to int
        pattern = re.compile('\d+')  # pattern for all decimal digits
        # print(pattern.findall(bias_str))

        # hold the result in an object
        bias_sec = pattern.findall(bias_str)

    # typecast the indices to ints
    xmin = int(bias_sec[0])
    xmax = int(bias_sec[1])
    ymin = int(bias_sec[2])
    ymax = int(bias_sec[3])

    # retrieve the slice of data corresponding to the bias section
    im_data = hdu.data
    bias_data = im_data[ymin:ymax, xmin:xmax]

    return bias_data


def sky_subtract(im_data, mask=None, mask_sources=True, **kwargs):
    """
    Returns a background (sky) subtracted copy of image data.

    Does so by calculating the mean, median and std of the background, and
    subtracting the mean from the input data. By default, this masks any sources
    and sigma clips.

    Parameters
    ----------
    im_data : ndarray
        requires an image data array
    mask: ndarray (bool), optional
        Boolean array the same size as im_data. Entries in im_data
        corresponding to True values in mask will be ignored in
        calculations.
    mask_sources: bool, optional
        If True, a source mask is automatically generated using photutils. This
        option allows this feature to be disabled, in case photutils is not
        available, or not working correctly.
    kwargs: dict, optional
        Keyword arguments to be passed to image_stats

    Returns
    -------
    output_im: ndarray
        The background subtracted data.
    stats: tuple
        The mean, median, and standard deviation of the background

    See Also
    --------
    background_stats: returns statistics on the background of image data
    """

    # calculate bias using mean
    # clipped stats are used, just in case
    print('----------------------')
    print('Background statistics:')
    mean, median, std = image_stats(im_data, mask=mask, mask_sources=mask_sources, **kwargs)

    # subtract the mean from the image data
    output_im = im_data - mean
    print('Background stats after subtraction:')
    image_stats(output_im, mask=mask, mask_sources=mask_sources, **kwargs)

    # package the statistics for return
    stats = mean, median, std
    # return generated values
    return output_im, stats


def bias_subtract(hdu, bias_sec=None):  # pass header data unit.  REMEBER, this is pass-by-reference
    """
    Returns the bias subtracted data from a header data unit.

    Takes a header data unit, and finds the mean of the bias section.
    It then subtracts that mean from a copy of the image data, and returns
    the result.

    Parameters
    ----------
    hdu : fits header data unit
        Image data stored in a fits file
    bias_sec: array-like, optional
        This must be a list, tuple, or array of four numbers, that define the
        bias section of the image. The numbers can be ints, floats, or strings
        that represent ints. If a float is given, it will be truncated
        towards zero.

    Returns
    -------
    output_im : numpy array
        the image data array after bias subtraction

    """

    # Store the data from the HDU argument
    im_data = hdu.data

    # check if parameters give the bias section, if not, automatically get it
    if bias_sec is None:
        # pull the bias section information from the header readout.
        bias_str = hdu.header['BIASSEC']
        print('Bias Section is ' + bias_str)
        # print(type(bias_str))
        # slice the string, for converting to int
        pattern = re.compile('\d+')  # pattern for all decimal digits
        print(pattern.findall(bias_str))

        # hold the result in an object
        bias_sec = pattern.findall(bias_str)

    # Bias section data
    # image is not indexed the same as python.
    # Image indexes (x,y), from lower left
    # python indexes (y,x)

    xmin = int(bias_sec[0])
    xmax = int(bias_sec[1])
    ymin = int(bias_sec[2])
    ymax = int(bias_sec[3])

    bias_data = im_data[ymin:ymax, xmin:xmax]

    # Calculate the bias, using clipped statistics in case of cosmic ray events, and print the 		#results
    bias_mean, bias_median, bias_std = sigma_clipped_stats(bias_data, sigma=3.0, iters=5)
    print('Bias mean: ' + str(bias_mean))
    print('Bias median: ' + str(bias_median))
    print('Bias standerd deviation: ' + str(bias_std))

    # calculate and print the bias area statistics, for reference.  DISABLED
    # print('Bias area after subtraction \n Mean: ')
    output_im = im_data - bias_mean

    return output_im


def display_data(imdata, **kwargs):
    """
    A wrapper for matplotlib.pyplot.imshow, for displaying image data.

    Parameters
    ----------
    imdata: ndarray
        image data array to be displayed
    """
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.figure()
    plt.imshow(imdata, norm=norm, origin='lower', cmap='viridis', **kwargs)
    plt.colorbar()
    plt.show()


def sigma_clipped_frame_average(filename_list, path='.', writeto_filename=None, sigma=3.0, iters=5, **kwargs):
    """
    This function calculates a sigma clipped average of images stored in fits
    files, by accessing them via a list of filenames.

    This function averages frame data that has been stored on disk in fits
    files. It does so by taking a list of file names. If a file path is not
    specified, defaults to the current working directory. The result can be
    written to file by specifying a filename, and is saved to the same
    directory as the source files. This function presumes that each fits file
    contains a Header Data Unit with extensions: the averages are computed on a
    per extension basis, and returned as a list of data arrays.

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
    path: str, optional
        A string that contains the file path that contains the filenames given.
        Defaults to the current working directory.
    writeto_filename: str or None, optional
        String that contains a file name to write the result to as HDU with
        extensions, in fits format. If None, does not write to file.
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
        hdul_list = [fits_stack.enter_context(fits.open(path + '/' + fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]
    # clean up the mess

    if writeto_filename is not None:
        frame_clipped_average_data, frame_clipped_average_error, frame_included_pixel_tracker = \
            _sigma_clipped_frame_average(image_stack, sigma=sigma, iters=iters,**kwargs)

        comment_string = 'Changed data to the sigma clipped average of ' + str(len(filename_list)) + \
                         ' zero frames, of which this is the first'

        _write_average_data_to_file(frame_clipped_average_data, writeto_filename, source_filename_list=filename_list,
                                    file_path=path, comment_string=comment_string)

        # print('Amount of garbage:')
        # print(gc.collect())

        return frame_clipped_average_data, frame_clipped_average_error, frame_included_pixel_tracker
    else:
        return _sigma_clipped_frame_average(image_stack, sigma=sigma, iters=iters,**kwargs)


def frame_average(filename_list, path='.', writeto_filename=None):
    """
    This function returns an average and error estimate of images stored in fits files, by
    accessing them via a list of filenames.

    This function averages frame data that has been stored on disk in fits
    files. It does so by taking a list of file names. If a file path is not
    specified, defaults to the current working directory. The result can be
    written to file by specifying a filename, and is saved to the same
    directory as the source files. This function presumes that each fits file
    contains a Header Data Unit with extensions. Data is taken from the
    extensions, stacked, and the averages computed. This produces a list of
    data arrays corresponding to the HDU data arrays.

    The error is calculated by taking the standard deviation of the data used
    to calculate the mean, on a pixel by pixel basis, and then dividing by the
    square root of the number of frames. This produces a list of arrays the
    same shape as the frame data.

    Parameters
    ----------
    filename_list: list
        List of strings containing the file names of the frames to be averaged.
    path: str, optional
        A string that contains the file path that contains the filenames given.
        Defaults to the current working directory.
    writeto_filename: str or None, optional
        String that contains a file name to write the result to as HDU with
        extensions, in fits format. If None, does not write to file.

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
        hdul_list = [fits_stack.enter_context(fits.open(path + '/' + fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    if writeto_filename is not None:
        frame_average_data, frame_average_error = _frame_mean(image_stack)
        comment_string = 'Changed data to the average of ' + str(len(filename_list)) + \
                         ' zero frames, of which this is the first'

        _write_average_data_to_file(frame_average_data, writeto_filename, source_filename_list=filename_list,
                                    file_path=path, comment_string=comment_string)
        return frame_average_data, frame_average_error
    else:
        return _frame_mean(image_stack)


def frame_median(filename_list, path='.', writeto_filename=None):
    """
    This function returns the median of several frames.

    The frames are accessed from disk, by passing the file names. If a file
    path is not specified, defaults to current working directory. The result
    can be written to a fits file by specifying a file name, and is saved to
    the same directory as the source files. This function presumes that each
    fits file contains a Header Data Unit with extensions. The data is stacked,
    and the median calculated on a per pixel basis. The median is intended to
    be a fast way of getting a bias frame, with some robustness against any
    outliers like cosmic rays.

    The error is calculated by finding the standard error on the mean, and
    multiplying by 1.253. Note that this only gives accurate error estimates
    for data that is normally distributed. This should be kept in mind when
    reviewing errors on data sets that contain outliers.

    Parameters
    ----------
    filename_list: list
        List of strings containing the file names of the frames to be averaged.
    path: str, optional
        A string that contains the file path for the file names given. Defaults
        to the current working directory.
    writeto_filename: str or None, optional
        String that contains a file name to write the result to as HDU with
        extensions, in fits format. If None, does not write to file.

    Returns
    -------
    frame_median: list
        A list of ndarrays containing image data after taking the median.
    pixel_error: list
        A list of ndarrays corresponding to frame_median, that contains the
        per-pixel error on the average.
    """

    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(path + '/' + fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    if writeto_filename is not None:
        frame_median_data, frame_median_error = _frame_median(image_stack)
        comment_string = 'Changed data to the median of' + str(len(filename_list)) + \
                         'zero frames, of which this is the first'

        _write_average_data_to_file(frame_median_data, writeto_filename, source_filename_list=filename_list,
                                    file_path=path, comment_string=comment_string)
        return frame_median_data, frame_median_error
    else:
        return _frame_median(image_stack)



def frame_subtract(minuend, subtrahend, display_in_ds9=False, write_to=None):
    """
    This function subtracts one HDUList from another, and returns the resultant
     data.

    This function expects a primary Header Data Unit and HDU extensions in a
    list. It can handle a single Primary HDU, but in that case expects a list
    of one. The HDULists can be passes directly as parameters using Astropy
    pyfits, as a filename string, or as an instance of pyds9.DS9. In the
    latter cases, the HDUlist will be opened from file or retrieved from DS9.
    It should be noted that DS9 will only pass a list of one HDU, the HDU
    that is loaded in the current frame.

    Parameters
    ----------
    minuend: HDUList, filename, or DS9 instance
        The source of the data to be subtracted from.
    subtrahend: HDUList, filename, or DS9 instance
        The source of the data to be subtracted.
    display_in_ds9: bool, optional
        If True, the result will be displayed in a Display instance of DS9,
        in a new frame. If the Display instance of DS9 is not already
        running, one will opened.
    write_to: str
        (To be implemented) Name of file to write result to. Will create a
        new file if one does not already exist.

    Returns
    -------
    difference:
        array containing the data after subtraction

    Examples
    --------
    >>> filename1 = '/home/lee/Documents/k4m_160531_050920_ori.fits.fz'
    >>> filename2 = '/home/lee/Documents/k4m_161228_132947_dri.fits.fz'
    >>> with fits.open(filename1) as minuend_hdul, fits.open(filename2) as subtrahend_hdul:
    ...     data_list = frame_subtract(minuend_hdul, subtrahend_hdul, display_in_ds9=False)
    >>> data_list[0]
    array([[4.700e+01, 5.700e+01, 1.049e+04, ..., 3.200e+01, 4.300e+01,
            3.800e+01],
           [3.400e+01, 2.800e+01, 3.225e+03, ..., 1.800e+01, 2.400e+01,
            3.900e+01],
           [3.300e+01, 1.800e+01, 2.774e+03, ..., 2.400e+01, 4.700e+01,
            3.900e+01],
           ...,
           [3.800e+01, 1.600e+01, 2.500e+01, ..., 2.100e+01, 3.300e+01,
            2.400e+01],
           [1.000e+01, 2.000e+01, 1.800e+01, ..., 1.500e+01, 1.300e+01,
            1.900e+01],
           [1.300e+01, 1.200e+01, 8.000e+00, ..., 9.000e+00, 1.800e+01,
            2.200e+01]])

    """
    if type(minuend) is pyds9.DS9:
        # if argument is a ds9 instance, open the current frame as an hdul
        minuend_hdul = minuend.get_pyfits()
    elif type(minuend) is str:
        # presume a string is a filename
        with fits.open(minuend) as hdul:
            # make a copy of the hdul from file
            minuend_hdul = fits.HDUList([hdu.copy() for hdu in hdul])
    else:
        # default cause, just pass the parameter directly
        minuend_hdul = minuend

    if type(subtrahend) is pyds9.DS9:
        # if argument is a ds9 instance, open the current frame as an hdul
        subtrahend_hdul = subtrahend.get_pyfits()
    elif type(subtrahend) is str:
        # presume a string is a filename
        with fits.open(subtrahend) as hdul:
            # make a copy of the hdul from file
            subtrahend_hdul = fits.HDUList([hdu.copy() for hdu in hdul])
    else:
        # default case, just pass directly
        subtrahend_hdul = subtrahend

    # special, but common use case:
    # one frame is open in DS9, and the other is on file
    # This requires matching the correct extension from file to the extension open in DS9
    # (this is mostly due to DS9 only passing primary HDUs)

    # handling the cases separately, because they require modifying different variables
    # first case: minuend is DS9.
    if type(minuend) is pyds9.DS9 and type(minuend) is not type(subtrahend):
        # retrieve extension name from DS9
        pattern = re.compile(r'(?<=\[).*(?=\])')  # pattern that retrieves everything between '[]'
        extension_name = pattern.search(minuend.get('file'))[0]

        subtrahend_hdul = [_find_hdu_extension(extension_name, subtrahend_hdul)]

    # second case: subtrahend is DS9. This case is not expected, but is included for robustness
    if type(subtrahend) is pyds9.DS9 and type(subtrahend) is not type(minuend):
        # retrieve extension name from DS9
        pattern = re.compile(r'(?<=\[).*(?=\])')  # pattern that retrieves everything between '[]'
        extension_name = pattern.search(subtrahend.get('file'))[0]

        # set minuend to be the corresponding extension
        minuend_hdul = [_find_hdu_extension(extension_name, minuend_hdul)]

    # calculate difference
    difference = _frame_subtract(minuend_hdul, subtrahend_hdul)

    if display_in_ds9:
        display = pyds9.DS9(target='Display', start='-title Display')
        for array in difference:
            display.set('frame new')
            display.set_np2arr(array)
    return difference


def get_filenames(path='.', extension=None, pattern=None, identifiers=None, include_path=False):
    """
    Retrieves a list containing the filenames from a target directory.

    By default this retrieves all entries in the directory, hereafter referred
    to as filenames. Specific file types or folders can be retrieved by
    filtering by extension, a portion of the filename, or a list of
    identifiers. Also supports Regular Expression patterns.

    Parameters
    ----------
    path: string or path-like object, optional
        The path from which to retrieve the filepaths. Default behavior is to
        list the current directory.
    extension: string, optional
        The extension of the filenames to be retrieved. This works by comparing
        the end of the entry names to the string specified, so it need not be
        an extension, merely the end of entry-name you wish to retrieve.
    pattern: string, optional
        A pattern by which to filter what filenames and entries are returned. This
        can be the whole filename, or just a portion. Alternately, a Regular
        Expression can be provided. All entries that match in this fashion will be
        returned.
    identifiers: list or tuple, optional
        Instead of a single string, a list of strings or numbers can be
        provided. The list can contain whole filenames, or just portions of the
        filenames. If a number is passed, the number is converted to a string.
        This does not support Regular Expressions.
    include_path: bool, optional
        If True, the filenames are returned with the path appended to them.

    Returns
    -------
    filename_list: a list of strings that contain the filenames after filtering

    See Also
    --------
    os.listdir: returns the names of all entries in a directory
    re: module that supports regular expression matching operations for python

    Examples
    --------

    Retrieve a list of file names using a file extension, and a Regex pattern.

    >>> biasframe_path = '/home/lee/Documents/bias_frames'
    >>> extension = '.fits.fz'
    >>> pattern = '(?=.*k4m)'  # look-ahead regex pattern that checks for 'k4m'
    >>> fits_list = get_filenames(biasframe_path, extension=extension, pattern=pattern)
    >>> print(fits_list)
    ['k4m_180211_231642_zri.fits.fz', 'k4m_180211_225927_zri.fits.fz', 'k4m_180211_231335_zri.fits.fz', 'k4m_180211_230119_zri.fits.fz', 'k4m_180211_231949_zri.fits.fz', 'k4m_180211_231834_zri.fits.fz', 'k4m_180211_231449_zri.fits.fz', 'k4m_180211_231604_zri.fits.fz', 'k4m_180211_231527_zri.fits.fz', 'k4m_180211_231911_zri.fits.fz', 'k4m_180211_230004_zri.fits.fz', 'k4m_180211_231719_zri.fits.fz', 'k4m_180211_231412_zri.fits.fz', 'k4m_180211_230042_zri.fits.fz', 'k4m_180211_231756_zri.fits.fz']

    Retrieve a list of file names using a sequence of numbers

    >>> identifiers = np.arange(190800, 190900)
    >>> fits_list = get_filenames(biasframe_path, identifiers=identifiers)
    >>> fits_list
    ['c4d_170331_190830_zri.fits.fz', 'c4d_170331_190853_zri.fits.fz']


    """
    # retrieve all filenames from the directory
    filename_list = os.listdir(path)

    # keep all filenames with the proper extension
    if extension is not None:
        filename_list = [filename for filename in filename_list if
                         filename[-len(extension):] == extension]

    # keep all filenames that match the pattern
    if pattern is not None:
        filename_list = [filename for filename in filename_list if re.search(pattern, filename)]

    # keep all filenames that match the identifiers provided
    if identifiers is not None:
        storage_list = []
        try:
            for ident in identifiers:
                storage_list.extend([filename for filename in filename_list if str(ident) in filename])

        except TypeError:
            print(identifiers, 'is not iterable')
        else:
            filename_list = storage_list

    if include_path:
        filename_list = [path + '/' + filename for filename in filename_list]

    return filename_list


