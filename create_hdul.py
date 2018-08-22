
"""
A function for generating a new hdul, after performing manipulations on exsisting hduls
"""

import os
import re

from astropy.io import fits


def get_filenames(path, extension='', pattern=''):
    # retrieve all filenames from the directory
    filename_list = os.listdir(path)

    # convert extension and pattern to raw strings
    r_extension = "%r"%extension
    # ensure all filenames have the proper extension

    fits_list = [path + '/' + filename for filename in filename_list if
                 re.search(extension+r'$', filename) and re.search(pattern, filename)]

    return fits_list

biasframe_path = '/home/lee/Documents/bias_frames'
extension = r'\.fits\.fz'
pattern = '(?=.*k4m)'

# retrieve the filenames
fits_list = get_filenames(biasframe_path, extension=extension, pattern=pattern)

# method for making a copy of a hdul
hdul = fits.open(fits_list[0])
hdu_list = [hdu.copy() for hdu in hdul]
hdulcopy = fits.HDUList(hdu_list)
hdul.close()

# method for combining two headers
hdul = fits.open(fits_list[1])
primary1 = hdulcopy[0]
primary2 = hdul[0]
# combine headers
primary1.header.extend(primary2.header)
print(primary1.header['DATE-OBS', 0])
print(primary1.header['DATE-OBS', 1])
print(primary2.header['DATE-OBS'])
hdul.close()
