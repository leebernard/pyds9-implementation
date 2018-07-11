

"""This file is for developing region statistical tools
"""
import numpy as np
from ccd_tools import *
# from ccd_tools import get_regions_from_ds9
from astropy.stats import sigma_clip, sigma_clipped_stats
# make sure to fix the implicit import


def get_single_region_from_ds9(get_data=True, ds9=None):
    """This function gets the first single valid box region selected in ds9

    """
    import pyds9
    import re

    # if a ds9 target is not specified, make one
    if not ds9:
        ds9 = pyds9.DS9()

    # set the region format to ds9 default, and coordinate system to image. This ensures the format is standardized.
    # image format is required to properly index the data array.
    ds9.set('regions format ds9')
    ds9.set('regions system image')

    # get selected regions info
    raw_string = ds9.get('regions selected')
    print(raw_string)

    # transform string into list that is organized by lines
    pattern = re.compile('.+')
    str_list = pattern.findall(raw_string)

    try:
        while not re.match('box', str_list[0]):
            if re.match('# tile', str_list[0]):
                print('')
            print(str_list.pop(0))
    except IndexError:
        message = 'No valid region selected in DS9. Please select a valid box region.'

def region_mean():
    """Calculates the mean of the selected regions in DS9

    Returns
    -------
    region_mean_list: List
        List of the mean values of the regions selected in DS9
    """

    region_data_list = [region.data for region in get_regions_from_ds9()]

    region_mean_list = [np.mean(data) for data in region_data_list]

    return region_mean_list


def region_median():
    """Calculates the median of the selected regions in DS9

    Returns
    -------
    region_median_list: list
        List of the median values of selected regions
    """
    region_data_list = [region.data for region in get_regions_from_ds9()]

    region_median_list = [np.median(data) for data in region_data_list]

    return region_median_list


def region_std():
    """Calculates the standard deviation of the selected regions  in DS9

    Returns
    -------
    region_std_list: list
        List of standard deviation values of the selected regions
    """
    region_data_list = [region.data for region in get_regions_from_ds9()]

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


def region_stats(sigma_clip_=False, mask_sources=False, mask=None):
    if mask_sources:
        pass# generate mask
    if sigma_clip:
        pass# return sigma clipped statistics, with mask = mask


