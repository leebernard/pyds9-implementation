

# a function for parsing the region info pulled from SAOImage DS9 by pyds9's access routines
# this function returns a list of region objects.
# Each object has the DS9 canonical definition of the region, the array indices of the region, and the region data
# for memory/runtime management concerns, the region data feature can be suppressed by setting the optional argument
# get_data=False. This prevents the function from accessing the data held in DS9, significantly decreasing the resource
# consumption.
def parse_regions(get_data=True):

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

    # This class is for convenient packaging of the region data
    class Region:
        pass

    # print meta data
    print('Region Coordinate system:')
    print(region_system)
    print('Selected Regions:')

    # parse the meta data string
    # pattern is all sequences of digits that are terminated by a period
    pattern = re.compile('\d+(?=\.)') #(?<!\.)

    # The list for holding the region data. This is returned
    regions = []

    for region_str in str_list:
        print(region_str)  # print the region currently being parsed
        if re.match('box', region_str):

            region_def = pattern.findall(region_str)

            # current instance of a region
            current_region = Region()

            # region format
            current_region.format = region_system

            # region definition: orgin is lower left, given as x and y coord, with a width and a height
            x_coord = int(region_def[0])
            y_coord = int(region_def[1])
            width = int(region_def[2])
            height = int(region_def[3])

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





