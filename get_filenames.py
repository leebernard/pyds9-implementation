

"""
Retrieves a list of filenames from a target directory.

Intended to retrive fits files for frame averaging and bias subtraction. Uses a pattern to determine what files to get,
and checks file extension
"""

import os
import re

biasframe_path = '/home/lee/Documents/bias_frames'
extension = r'\.fits\.fz'
pattern = '(?=.*k4m)(?=.*231642)'

def get_filenames(path, extension='', pattern=''):
    # retrieve all filenames from the directory
    filename_list = os.listdir(biasframe_path)

    # convert extension and pattern to raw strings
    r_extension = "%r"%extension
    # ensure all filenames have the proper extension

    fits_list = [biasframe_path + '/' + filename for filename in filename_list if
                 re.search(extension+r'$', filename) and re.search(pattern, filename)]

    return fits_list

