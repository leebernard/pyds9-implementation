

# import needed packages
import pyds9
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
import warnings
import numpy as np
# from astropy.io import fits

import re
from astropy.stats import sigma_clipped_stats


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

    def stats(self, sigma_clip=False, mask=None, **kwargs):
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
        kwargs: dict
            Keyword argument
        Returns
        -------

        """
        if sigma_clip:
            mean, median, std = sigma_clipped_stats(self.data, mask, **kwargs)
        else:
            masked_data = np.ma.array(self.data, mask=mask)
            mean, median, std = np.ma.mean(masked_data), np.ma.median(masked_data), np.ma.std(masked_data)

        print('-----------------------')
        print('Region Data Statistics:')
        print(f'Mean: {mean:.2f}')
        print(f'Median: {median:.2f}')
        print(f'Std: {std:.2f}')
        # if there are stats on the bias and sky background, print those also
        if self.bias_sub:
            if sigma_clip:
                mean, median, std = sigma_clipped_stats(self.bias_sub, mask, **kwargs)
            else:
                masked_data = np.ma.array(self.bias_sub, mask=mask)
                mean, median, std = np.ma.mean(masked_data), np.ma.median(masked_data), np.ma.std(masked_data)
            print('---------------------')
            print('Bias Subtracted Data Statistics:')
            print(f'Mean: {mean:.2f}')
            print(f'Median: {median:.2f}')
            print(f'Std: {std:.2f}')
        if self.sky_sub:
            if sigma_clip:
                mean, median, std = sigma_clipped_stats(self.bias_sub, mask, **kwargs)
            else:
                masked_data = np.ma.array(self.bias_sub, mask=mask)
                mean, median, std = np.ma.mean(masked_data), np.ma.median(masked_data), np.ma.std(masked_data)
            print('---------------------')
            print('Sky Subtracted Data Statistics:')
            print(f'Mean: {self.sky_stats[0]:.2f}')
            print(f'Median: {self.sky_stats[1]:.2f}')
            print(f'Std: {self.sky_stats[2]:.2f}')

        # return mean, median, std

    def sky_subtract(self, mask=None, **kwargs):
        """
        This function calculates the background subtracted data, and stores it in
        the sky_sub attribute.

        Wrapper for sky_subtract. If bias subtracted data is available, it
        uses that. Otherwise, it uses data in the data attribute. If neither are
        available, raises an exception.

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
        if self.bias_sub is not None:
            self.sky_sub, self.sky_stats = sky_subtract(self.bias_sub, mask=mask, **kwargs)
        # if bias subtracted data is not found, used raw data
        elif self.data is not None:
            warnings.warn('Bias subtracted data not found. Using raw data.', category=UserWarning)
            self.sky_sub, self.sky_stats = sky_subtract(self.data, mask=mask, **kwargs)
        # if no data is available, throw an exception
        elif self.data is None:
            raise ValueError('Data attributes are unassigned')


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
    pattern = re.compile('.+')
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
    # pattern is all sequences of digits that may or may not contain a period
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
        plt.imshow(mask, origin='lower', cmap='viridis')
    return mask


def image_stats(indata, mask=None, mask_sources=False, sigma_clip=True, **kwargs):
    """
    Calculates statistics on the background of the image data given.

    Parameters
    ----------
    indata: array-like
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
        obj_mask = make_source_mask(indata, display_mask=True)
        # if both a mask is specified and object masking is called for, combine the masks
        if mask_sources and mask:
            mask = mask + obj_mask
        # otherwise, mask is just the object mask
        else:
            mask = obj_mask

    # calculate the stats, with optional masks and optional sigma clipping
    if sigma_clip:
        mean, median, std = sigma_clipped_stats(indata, mask=mask, **kwargs)
    else:
        # mask the data. If no mask is specified, mask is None
        masked_data = np.ma.array(indata, mask=mask)
        mean, median, std = np.ma.mean(masked_data), np.ma.median(masked_data), np.ma.std(masked_data)

    # display and return the values
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

    # subtract the bias from the image data
    output_im = im_data - mean
    print('Background stats after subtraction:')
    image_stats(output_im, mask=mask, mask_sources=mask_sources, **kwargs)

    # package the statistics for return
    stats = mean, median, std
    # return generated values
    return output_im, stats


def bias_subtract(HDU, bias_sec=None):  # pass header data unit.  REMEBER, this is pass-by-reference
    """
    Returns the bias subtracted data from a header data unit.

    Takes a header data unit, and finds the mean of the bias section.
    It then subtracts that mean from a copy of the image data, and returns
    the result.

    Parameters
    ----------
    HDU : fits header data unit
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
    im_data = HDU.data

    # check if parameters give the bias section, if not, automatically get it
    if bias_sec is None:
        # pull the bias section information from the header readout.
        bias_str = HDU.header['BIASSEC']
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


def get_image_data(ds9=None, datasec=None):

    # if a ds9 target is not specified, make one
    if ds9 is None:
        try:
            ds9 = pyds9.DS9()
        except ValueError:
            print('Specify a target DS9() instance to retrieve the bias from. '
                  'e.g:\n  d = pyds9.DS9(\'7f000101:43123\')\n  r = get_ds9_region(ds9=d)')
            raise
    hdulist = ds9.get_pyfits()

    if datasec is None:
        # determine what the valid data region is
        pattern = re.compile('\d+')  # pattern for all decimal digits
        datasec = pattern.findall(hdulist[0].header['DATASEC'])
        print('Data section: ', datasec)

    image_data = hdulist[0].data[int(datasec[0]):int(datasec[1]), int(datasec[2]):int(datasec[3])]

    return image_data

def get_multiple_ds9_regions(get_data=True, ds9=None):
    """a function for importing the region info from SAOImage DS9 by pyds9's access routines.

    Each object has the DS9 canonical definition of the region, the array indices of the region, and the region data
    for memory/runtime management concerns, the region data feature can be suppressed by setting the optional argument
    get_data=False. This prevents the function from accessing the data held in DS9, significantly decreasing the
    resource consumption.

    Parameters
    ----------
    get_data: bool, optional
        sets whether or not to include the data of the region
    ds9: DS9 object, optional
        Gives the function a ds9 object to pull the regions from. If none is indicated, creates an DS9 object using
        default methods.

    Returns
    -------
    regions: list
        list of region objects that are selected in SAOImage DS9
    """
    # pulls all regions into a list. 1st entry on the list is the frame name
    import pyds9
    import re
    # import numpy as np

    # if ds9 is not specified, create a ds9 object
    if not ds9:
        ds9 = pyds9.DS9()

    # set the region format to ds9 default, and coordinate system to image. This ensures the format is standardized.
    # image format is required to properly index the data array.
    ds9.set('regions format ds9')
    ds9.set('regions system image')

    # get selected regions info
    raw_string = ds9.get('regions selected')
    # print(raw_string)

    # transform string into list that is organized by lines
    pattern = re.compile('.+')

    str_list = pattern.findall(raw_string)

    # remove meta-meta data, first two entries
    del str_list[0:2]

    # failure condition: no regions selected
    try:
        # yank format
        region_system = str_list.pop(0)
    except IndexError:
        print('No region selected in DS9. Please select a region')
        return 1

    # retrieve frame data
    if get_data:
        frame_data = ds9.get_arr2np()

    # frame_name = 'current frame'

    # print meta data
    print('Region Coordinate system:')
    print(region_system)
    print('Selected Regions:')

    # parse the meta data string
    # pattern is all sequences of digits that may or may not contain a period
    pattern = re.compile('\d+\.?\d*')

    # The list for holding the region data. This is returned
    regions = []

    for region_str in str_list:
        print(region_str)  # print the region currently being parsed
        if re.match('box', region_str):

            region_def = pattern.findall(region_str)

            # current instance of a region
            current_region = Region()

            # save region definition
            current_region.region_def = region_str

            # region system
            current_region.system = region_system

            # region definition: origin is lower left, given as x and y coord, with a width and a height
            x_coord = float(region_def[0])
            y_coord = float(region_def[1])
            width = float(region_def[2])
            height = float(region_def[3])

            # region slicing data
            xmin = int(x_coord - width/2)
            xmax = int(x_coord + width/2)

            ymin = int(y_coord - height/2)
            ymax = int(y_coord + height/2)

            # retrieve region data by slicing the frame data array, and store it in the current region
            # This is determined by option get_data.
            if get_data:
                current_region.data = frame_data[ymin:ymax, xmin:xmax]

            # package all the meta data
            current_region.x_coord = x_coord
            current_region.y_coord = y_coord
            current_region.width = width
            current_region.height = height

            current_region.xmin = xmin
            current_region.xmax = xmax
            current_region.ymin = ymin
            current_region.ymax = ymax

            # store current region in the return list
            regions.append(current_region)

            print('Region resolved')

        else:
            print('Region is not a box!')  # error condition
    return regions
