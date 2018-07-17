

import pyds9
import ccd_tools


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