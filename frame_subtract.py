"""
A function for subtracting (and adding) frames from fits files
"""
import numpy as np
import pyds9
from astropy.io import fits
import re

# import timing


def _find_hdu_extension(extname, hdulist):
    # finds the first match in the hdulist for the extension name
    # This should work for any format of hdu list
    for hdu in hdulist:
        try:
            if re.match(extname, hdu.header['extname'], flags=re.IGNORECASE):
                hdu_match = hdu
                # once the extension is found, break
                break
        except KeyError:
            # catch the exception when no extension name exists, and keep moving
            print('No extension name.')
    else:
        # corresponding extension not found, raise error
        raise RuntimeError('A matching extension was not found.')

    return hdu_match


def _frame_subtract(minuend_hdul, subtrahend_hdul):
    return [minuend.data.astype('float64', casting='safe') - subtrahend.data.astype('float64', casting='safe')
            for minuend, subtrahend in zip(minuend_hdul, subtrahend_hdul)
            if minuend.data is not None and subtrahend.data is not None]


def frame_subtract(minuend, subtrahend, display=True, write_to=None):
    """
    This function subtracts one HDUList from another, and returns the resultant
     data.

    This function expects a primary Header Data Unit and HDU extensions in a
    list. It can handle a single Primary HDU, but in that case expects a list
    of one. The HDULists can be passes directly as parameters using Astropy
    pyfits, as a filename string, or as an instance of pyds9.DS9. In the
    latter cases, the HDUlist will be opened from file or retrieved from DS9.
    It should be noted that DS9 will only pass a list of one HDU, the HDU
    that is loaded in the current frame.

    Parameters
    ----------
    minuend: HDUList, filename, or DS9 instance
        The source of the data to be subtracted from.
    subtrahend: HDUList, filename, or DS9 instance
        The source of the data to be subtracted.
    display: bool, optional
        If True, the result will be displayed in a Display instance of DS9,
        in a new frame. If the Display instance of DS9 is not already
        running, one will opened.
    write_to: str
        (To be implemented) Name of file to write result to. Will create a
        new file if one does not already exist.

    Returns
    -------
    difference:
        array containing the data after subtraction
    """
    if type(minuend) is pyds9.DS9:
        # if argument is a ds9 instance, open the current frame as an hdul
        minuend_hdul = minuend.get_pyfits()
    elif type(minuend) is str:
        # presume a string is a filename
        with fits.open(minuend) as hdul:
            # make a copy of the hdul from file
            minuend_hdul = fits.HDUList([hdu.copy() for hdu in hdul])
    else:
        # default cause, just pass the parameter directly
        minuend_hdul = minuend

    if type(subtrahend) is pyds9.DS9:
        # if argument is a ds9 instance, open the current frame as an hdul
        subtrahend_hdul = subtrahend.get_pyfits()
    elif type(subtrahend) is str:
        # presume a string is a filename
        with fits.open(subtrahend) as hdul:
            # make a copy of the hdul from file
            subtrahend_hdul = fits.HDUList([hdu.copy() for hdu in hdul])
    else:
        # default case, just pass directly
        subtrahend_hdul = subtrahend

    # special, but common use case:
    # one frame is open in DS9, and the other is on file
    # This requires matching the correct extension from file to the extension open in DS9
    # (this is mostly due to DS9 only passing primary HDUs)

    # handling the cases separately, because they require modifying different variables
    # first case: minuend is DS9.
    if type(minuend) is pyds9.DS9 and type(minuend) is not type(subtrahend):
        # retrieve extension name from DS9
        pattern = re.compile(r'(?<=\[).*(?=\])')
        extension_name = pattern.search(minuend.get('file'))[0]

        subtrahend_hdul = [_find_hdu_extension(extension_name, subtrahend_hdul)]

    # second case: subtrahend is DS9. This case is not expected, but is included for robustness
    if type(subtrahend) is pyds9.DS9 and type(subtrahend) is not type(minuend):
        # retrieve extension name from DS9
        pattern = re.compile(r'(?<=\[).*(?=\])')
        extension_name = pattern.search(subtrahend.get('file'))[0]

        # set minuend to be the corresponding extension
        minuend_hdul = [_find_hdu_extension(extension_name, minuend_hdul)]

    # calculate difference
    difference = _frame_subtract(minuend_hdul, subtrahend_hdul)

    if display:
        display = pyds9.DS9(target='Display', start='-title Display')
        for array in difference:
            display.set('frame new')
            display.set_np2arr(array)
    return difference


# regex pattern for finding everything between brackets '[]'
# filename1 = '/home/lee/Documents/bias_frames/c4d_170331_202028_zri.fits.fz'
# filename2 = '/home/lee/Documents/bias_frames/c4d_170331_191113_zri.fits.fz'

filename1 = '/home/lee/Documents/k4m_160531_050920_ori.fits.fz'
filename2 = '/home/lee/Documents/k4m_161228_132947_dri.fits.fz'
biasframe_file = '/home/lee/Documents/bias_frames/testaverage.fits.fz'
# Open a new ds9 instance, or if already open, access it
display = pyds9.DS9(target='display', start='-title display')
ds9 = pyds9.DS9()
ds9.set('file ' + filename1)
with fits.open(filename1) as minuend_hdul, fits.open(filename2) as subtrahend_hdul:
    # minuend_hdul.info()
    # subtrahend_hdul.info()

    data_list = frame_subtract(minuend_hdul, subtrahend_hdul)

    difference = frame_subtract(ds9, filename2)
    # for data in data_list:
    #     display.set('frame new')
    #     display.set_np2arr(data)


