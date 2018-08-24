

"""
Retrieves a list of filenames from a target directory.

Intended to retrive fits files for frame averaging and bias subtraction. Uses a pattern to determine what files to get,
and checks file extension
"""

import os
import re

# path = '/home/lee/Documents/bias_frames'
# extension = '.fits.fz'
# pattern = '(?=.*k4m)'
# identifiers = [120]

def get_filenames(path, extension=None, pattern=None, identifiers=None, include_path=False):
    """
    Retrieves a list containing the filenames from a target directory.

    By default this retrieves all entries in the directory, hereafter referred
    to as filenames. Specific file types or folders can be retrieved by
    filtering by extension, a portion of the filename, or a list of
    identifiers. Also supports Regular Expression patterns.

    Parameters
    ----------
    path: string or path-like object
        The path from which to retrieve the filepaths
    extension: string, optional
        The extension of the filenames to be retrieved. This works by comparing
         the end of the entry names to the string specified, so it need not be
         an extension, merely the end of entry-name you wish to retrieve
    pattern: string, optional
        A pattern by which to filter what filenames and entries are returned. This
        can be the whole filename, or just a portion. Alternately, a Regular
        Expression can be provided. All entries that match in this fashion will be
        returned.
    identifiers: list or tuple, optional
        Instead of a single string, a list of strings or numbers can be
        provided. The list can contain whole filenames, or just portions of the
        filenames. This does not support Regular Expressions, but does support
        integers.
    include_path: bool, optional
        If True, the filenames are returned with the path appended to them.

    Returns
    -------
    filename_list: a list of strings that contain the filenames after filtering

    See Also
    --------
    os.listdir: returns the names of all entries in a directory
    re: module that supports regular expression matching operations for python
    """
    # retrieve all filenames from the directory
    filename_list = os.listdir(path)

    # keep all filenames with the proper extension
    if extension is not None:
        filename_list = [filename for filename in filename_list if
                         filename[-len(extension):] == extension]

    # keep all filenames that match the pattern
    if pattern is not None:
        filename_list = [filename for filename in filename_list if
                     re.search(pattern, filename)]

    # keep all filenames that match the identifiers provided
    if identifiers is not None:
        storage_list = []
        try:
            for ident in identifiers:
                storage_list.extend([filename for filename in filename_list if str(ident) in filename])

        except TypeError:
            print(identifiers, 'is not iterable')
        else:
            filename_list = storage_list

    if include_path:
        filename_list = [path + '/' + filename for filename in filename_list]

    return filename_list

