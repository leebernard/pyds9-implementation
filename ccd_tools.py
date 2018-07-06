


def bias_subtract(HDU, bias_sec=None):  # pass header data unit.  REMEBER, this is pass-by-reference
    """Takes a header data unit, find the bias data from BIASSEC, and performs bias calculations and subtraction.

    Parameters
    ----------
    HDU : fits header data unit
        Image data stored in a fits file
    bias_sec : tuple int, optional
        defines the area of the frame to be used to calculate the bias. If not specified, determines the bias from the
        header definition
    Returns
    -------
    output_im : numpy array
        the image data array after bias subtraction

    """
    # import needed packages
    # import numpy as np
    # from astropy.io import fits
    import re
    from astropy.stats import sigma_clipped_stats

    # Store the data from the HDU argument
    im_data = HDU.data

    # check if parameters give the bias section, if not, automatically get it
    if not bias_sec:
        # pull the bias section information from the header readout.
        Bias_Sec = HDU.header['BIASSEC']
        print('Bias Section is ' + Bias_Sec)
        # print(type(Bias_Sec))
        # slice the string, for converting to int
        pattern = re.compile('\d+')  # pattern for all decimal digits
        print(pattern.findall(Bias_Sec))

        # hold the result in an object
        bias_sec = pattern.findall(Bias_Sec)

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



def background_subtract(im_data):
    """calculates background using a mask routine from photutils.

    Parameters
    ----------
    im_data : numpy array
        requires an image data array

    returns
    ------
    output_im: numpy array
        The background subtracted data.
    mask: numpy array bool
        The mask used to shield the object
    """
    # import numpy as np
    # from astropy.io import fits

    # store the data from the HDU argument
    # im_data = HDU.data

    # Generate mask
    from photutils import make_source_mask
    from astropy.stats import sigma_clipped_stats
    mask = make_source_mask(im_data, snr=2, npixels=5, dilate_size=11)

    # calculate bias using mean
    # clipped stats are used, just in case
    mean, median, std = sigma_clipped_stats(im_data, sigma=3.0, mask=mask)
    print('Background mean: ' + str(mean))
    print('Background median: ' + str(median))
    print('Background standerd deviation: ' + str(std))

    output_im = im_data - mean

    return output_im, mask, std



def get_regions(get_data=True):
    """a function for importing the region info from SAOImage DS9 by pyds9's access routines.

    Each object has the DS9 canonical definition of the region, the array indices of the region, and the region data
    for memory/runtime management concerns, the region data feature can be suppressed by setting the optional argument
    get_data=False. This prevents the function from accessing the data held in DS9, significantly decreasing the
    resource consumption.

    Parameters
    ----------
    get_data: bool, optional
        sets whether or not to include the data of the region

    Returns
    -------
    regions: list
        list of region objects that are selected in SAOImage DS9
    """
    # pulls all regions into a list. 1st entry on the list is the frame name
    import pyds9
    import re
    # import numpy as np


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


    class Region:
        """
        This class is for convenient packaging of the region data.

        quick guide
        -----------
        self.xmin, self.ymin:
            gives array location of lower left corner of region
        self.xmax, self.ymax:
            gives array location of upper right corner of region
        self.data:
            data array selected by region

        Attributes
        ----------
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
        data: int, optional
            data array of the defined region, pulled directly from SAOImage DS9. This can be disabled by the get_data
            flag
        """
        pass

    # print meta data
    print('Region Coordinate system:')
    print(region_system)
    print('Selected Regions:')

    # parse the meta data string
    # pattern is all sequences of digits that are terminated by a period
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

            # region definition: orgin is lower left, given as x and y coord, with a width and a height
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

