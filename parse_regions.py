# a function for parsing the region info pulled from DS9 by pyds9's access routines

# pull a sample region
import pyds9
import re


ds9 = pyds9.DS9()

# set format to image. This ensures the format is standardized
ds9.set('regions system image')

# get selected regions info
raw_string = ds9.get('regions selected')
# print(raw_string)

# re pattern that pulls all strings that do not include an newline char
pattern = re.compile('.+')

match = pattern.findall(raw_string)


# remove meta-meta data, first two entries
del match[0:2]


print('Meta data format')
print(match.pop(0))
for str in match:
     print(str)
