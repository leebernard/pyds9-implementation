# a function for parsing the region info pulled from DS9 by pyds9's access routines

# pull a sample region
import pyds9
import re
import numpy as np

ds9 = pyds9.DS9()

# set format to image. This ensures the format is standardized
ds9.set('regions system image')

# get selected regions info
raw_string = ds9.get('regions selected')
# print(raw_string)

# transform string into list that is organized by lines
pattern = re.compile('.+')

str_list = pattern.findall(raw_string)


# remove meta-meta data, first two entries
del str_list[0:2]

# yank format
region_format = str_list.pop(0)

# print meta data
print('Meta data format')
print(region_format)
print('Selected Regions')
for str in str_list:
     print(str)


frame_data = ds9.get_arr2np()
# parse the meta data string
n=0
if re.match('box',str_list[n]):
    print('Region resolved')
    pattern = re.compile('\d+')

    region_def = pattern.findall(str_list[0])

    # region definition: orgin is lower left, given as x and y coord, with a width and a height
    x_coord = int(region_def[0])
    y_coord = int(region_def[1])
    width = int(region_def[2])
    height = int(region_def[3])

    xmin = x_coord
    xmax = xmin + width

    ymin = y_coord
    ymax = ymin + height

    # pull current region
    current_region = frame_data[ymin:ymax, xmin:xmax]


else:
    print('Region is not a box!')  # error condition
# retrieve frame data



