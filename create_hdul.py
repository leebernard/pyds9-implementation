
"""
A function for generating a new hdul, after performing manipulations on exsisting hduls
"""

import os
import re

from astropy.io import fits
from get_filenames import get_filenames


biasframe_path = '/home/lee/Documents/bias_frames'
extension = '.fits.fz'
pattern = '(?=.*k4m)'
identifiers = ['231642', '225927', '231335']
# retrieve the filenames
fits_list = get_filenames(biasframe_path, extension=extension, pattern=pattern, include_path=True)

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
