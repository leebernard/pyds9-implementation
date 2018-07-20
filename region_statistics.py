

"""
This file is for developing region statistical tools
"""
import numpy as np
import pyds9
import re
from ccd_tools import *
# from ccd_tools import get_multiple_ds9_regions
from astropy.stats import sigma_clip, sigma_clipped_stats
# make sure to fix the implicit import


def get_ds9_region(get_data=True, ds9=None, tiled=False):
    """
    This function gets the first single valid box region selected in ds9

    Parameters
    ----------
    get_data: bool, optional
        If True, does not retrieve data from DS9. This to reduce the resource
        requirements if data is being handled separately.
    ds9: DS9 object, optional
        optional parameter to specify a DS9 target

    Returns
    -------
    region: Region object

    """


    # check if ds9 is accesible
    if pyds9.ds9_targets() is None:
        input('DS9 target not found. Please start/restart DS9, then press enter')

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

    if get_data:
        with ds9.get_arr2np() as frame_data:

            region.data = frame_data[ymin:ymax, xmin:xmax]


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
    region.region_def = raw_string

    return region


def region_stats(self, sigma_clip_=False, **kwargs):

    if sigma_clip_:
        return




def multiple_region_mean():
    """Calculates the mean of the selected regions in DS9

    Returns
    -------
    region_mean_list: List
        List of the mean values of the regions selected in DS9
    """

    region_data_list = [region.data for region in get_multiple_ds9_regions()]

    region_mean_list = [np.mean(data) for data in region_data_list]

    return region_mean_list


def multiple_region_median():
    """Calculates the median of the selected regions in DS9

    Returns
    -------
    region_median_list: list
        List of the median values of selected regions
    """
    region_data_list = [region.data for region in get_multiple_ds9_regions()]

    region_median_list = [np.median(data) for data in region_data_list]

    return region_median_list


def multiple_region_std():
    """Calculates the standard deviation of the selected regions  in DS9

    Returns
    -------
    region_std_list: list
        List of standard deviation values of the selected regions
    """
    region_data_list = [region.data for region in get_multiple_ds9_regions()]

    region_std_list = [np.std(data) for data in region_data_list]

    return region_std_list


def sort_regions(region_list):
    """Sorts a list of region objects by their distance from the lower left corner in the array

    sorts by the center of the region
    
    Parameters
    ----------
    region_list: list
        List of Region objects to be sorted
    """
    region_list.sort(key=lambda region: np.sqrt(region.x_coord**2 + region.y_coord**2))


def print_region_def(region_list):
    [print(region.region_def) for region in region_list]


