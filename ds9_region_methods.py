
"""
Development space
"""

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

def get_ds9_region(get_data=True, ds9=None, tiled=False):
    """This function gets the first single valid box region selected in ds9

    Parameters
    ----------
    get_data: bool, optional
        flag to disable the data retrieval from DS9
    ds9: DS9 object, optional
        optional parameter to specify a DS9 target

    Returns
    -------
    current_region: Region object


    """

    import re

    # check if ds9 is accesible
    if pyds9.ds9_targets() is None:
        input('DS9 target not found. Please start/restart DS9, then press enter')

    # if a ds9 target is not specified, make one
    if ds9 is None:
        ds9 = pyds9.DS9()

    # set the region format to ds9 default, and coordinate system to image. This ensures the format is standardized.
    # image format is required to properly index the data array.
    ds9.set('regions format ds9')
    ds9.set('regions system image')

    # get selected region info
    raw_string = ds9.get('regions selected')
    print(raw_string)

    # transform string into list that is organized by lines
    pattern = re.compile('.+')
    str_list = pattern.findall(raw_string)

    try:
        while not re.match('box', str_list[0]):
            if re.match('# tile', str_list[0]):
                print('Tile mode detected')
                tiled = True

            print(str_list.pop(0))
    except IndexError:
        message = 'No valid region found. Please select a valid box region.'
        print(message)

    # parse the meta data string
    # pattern is all sequences of digits that may or may not contain a period
    pattern = re.compile('\d+\.?\d*')
    region_def = pattern.findall(str_list[0])  # should probably add an exception test here


    # make a region object to hold all the data
    current_region = ccd_tools.Region()

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

    if get_data:
        with ds9.get_arr2np() as frame_data:
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

    return current_region



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
