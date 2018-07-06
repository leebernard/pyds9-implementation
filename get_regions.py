

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
        raise


    # retrieve frame data
    if get_data:
        frame_data = ds9.get_arr2np()

    # frame_name = 'current frame'

    class _Region:
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
            current_region = _Region()

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

