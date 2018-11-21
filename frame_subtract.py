"""
A function for subtracting (and adding) frames from fits files
"""
import numpy as np
import pyds9
import datetime
import re
import warnings

from astropy.io import fits


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


def _copy_hdul(hdul):
    # iterate through the HDUList, copying each header data unit
    hdu_list = [hdu.copy() for hdu in hdul]
    # generate and return a new HDUList
    return fits.HDUList(hdu_list)


def _write_difference_to_file(data_list, writeto_filename, minuend_hdul, file_path='.', comment_string=None):
    """
    This is a specific method for frame subtraction.


    Parameters
    ----------
    data_list
    writeto_filename
    source_file_list
    file_path
    comment_string

    """

    # add a None to the start of the data list, for the primary HDU
    data_list.insert(0, None)
    # make generator for modifying the fits file
    hdul_generator = (hdu for hdu in minuend_hdul)

    for hdu, diff_data in zip(hdul_generator, data_list):
        hdu.data = diff_data
        hdu.header.set('OBSTYPE', 'subtracted')
        # input('press Enter to continue')
        if type(hdu) is fits.hdu.image.PrimaryHDU:
            hdu.header.add_comment(str(datetime.date.today()))
            hdu.header.add_comment('Modified OBSTYPE')
            if comment_string:
                hdu.header.add_comment(comment_string)

    minuend_hdul.writeto(file_path + '/' + writeto_filename)


def frame_subtract(minuend, subtrahend, file_path='.', display_in_ds9=False, write_to=None):
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
    file_path: string, optional
        The directory that contains the file names. Default is the current directory.
    display_in_ds9: bool, optional
        If True, the result will be displayed in a Display instance of DS9,
        in a new frame. If the Display instance of DS9 is not already
        running, one will opened.
    write_to: str
        Name of file to write result to. It will create a new file if one does not
        already exist. Due to inconsistencies between how python and SAOImage DS9 handles
        header data units, this is disallowed when DS9 is one or both of the sources.

    Returns
    -------
    difference:
        array containing the data after subtraction

    Examples
    --------
    >>> filename1 = '/home/lee/Documents/k4m_160531_050920_ori.fits.fz'
    >>> filename2 = '/home/lee/Documents/k4m_161228_132947_dri.fits.fz'
    >>> with fits.open(filename1) as minuend_hdul, fits.open(filename2) as subtrahend_hdul:
    ...     data_list = frame_subtract(minuend_hdul, subtrahend_hdul, display_in_ds9=False)
    >>> data_list[0]
    array([[4.700e+01, 5.700e+01, 1.049e+04, ..., 3.200e+01, 4.300e+01,
            3.800e+01],
           [3.400e+01, 2.800e+01, 3.225e+03, ..., 1.800e+01, 2.400e+01,
            3.900e+01],
           [3.300e+01, 1.800e+01, 2.774e+03, ..., 2.400e+01, 4.700e+01,
            3.900e+01],
           ...,
           [3.800e+01, 1.600e+01, 2.500e+01, ..., 2.100e+01, 3.300e+01,
            2.400e+01],
           [1.000e+01, 2.000e+01, 1.800e+01, ..., 1.500e+01, 1.300e+01,
            1.900e+01],
           [1.300e+01, 1.200e+01, 8.000e+00, ..., 9.000e+00, 1.800e+01,
            2.200e+01]])

    >>> source_directory = '/home/lee/Documents'
    ... filename1 = 'k4m_160531_050920_ori.fits.fz'
    ... filename2 = 'k4m_161228_132947_dri.fits.fz'
    >>> __ = frame_subtract(filename1, filename2, file_path=source_directory,
    ...                     display_in_ds9=True, write_to='test_result.fits.fz')
    >>> test_result = fits.open('/home/lee/Documents/test_result.fits.fz')
    >>> test_result
    [<astropy.io.fits.hdu.image.PrimaryHDU object at 0x7fb4c58dab38>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49a914630>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49e90f320>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49a9a8be0>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49aa35f28>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb494611208>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49a91a400>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb496bcaef0>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49494d668>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49490f7f0>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb4948c6390>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb494888c88>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb4948487b8>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb49480bc18>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb4947c44a8>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb498433908>,
    <astropy.io.fits.hdu.compressed.CompImageHDU object at 0x7fb498199d30>]
    >>> test_result.close()
    """
    if type(minuend) is pyds9.DS9:
        # if argument is a ds9 instance, open the current frame as an hdul
        minuend_hdul = minuend.get_pyfits()
    elif type(minuend) is str:
        # presume a string is a filename
        # open the filename in the indicated directory
        with fits.open(file_path+'/'+minuend) as hdul:
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
        # open the filename in the indicated directory
        with fits.open(file_path+'/'+subtrahend) as hdul:
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
        pattern = re.compile(r'(?<=\[).*(?=\])')  # pattern that retrieves everything between '[]'
        extension_name = pattern.search(minuend.get('file'))[0]

        subtrahend_hdul = [_find_hdu_extension(extension_name, subtrahend_hdul)]

    # second case: subtrahend is DS9. This case is not expected, but is included for robustness
    if type(subtrahend) is pyds9.DS9 and type(subtrahend) is not type(minuend):
        # retrieve extension name from DS9
        pattern = re.compile(r'(?<=\[).*(?=\])')  # pattern that retrieves everything between '[]'
        extension_name = pattern.search(subtrahend.get('file'))[0]

        # set minuend to be the corresponding extension
        minuend_hdul = [_find_hdu_extension(extension_name, minuend_hdul)]

    # calculate difference
    difference = _frame_subtract(minuend_hdul, subtrahend_hdul)

    if display_in_ds9:
        display = pyds9.DS9(target='Display', start='-title Display')
        for array in difference:
            display.set('frame new')
            display.set_np2arr(array)

    if write_to:
        # if an argument is DS9, throw an expection
        if type(minuend) == pyds9.DS9 or type(subtrahend) == pyds9.DS9:
            raise TypeError('Saving results automatically is disallowed when a source is DS9.')

        # if one or both arguments are HDUs, throw a warning, and continue
        if type(minuend) == fits.hdu.hdulist.HDUList:
            warnings.warn('Source file name for the minuend is unknown.', category=UserWarning)
            minuend_source = 'unknown'
        else:
            minuend_source = minuend

        if type(subtrahend) == fits.hdu.hdulist.HDUList:
            warnings.warn('Source file name for the subtrahend is unknown', category=UserWarning)
            subtrahend_source = 'unknown'
        else:
            subtrahend_source = subtrahend

        # generate a comment string that updates the header with the source file names
        comment_string = 'Result of subtraction of '+subtrahend_source+' from '+minuend_source+'.'

        # save the result to file
        _write_difference_to_file(difference, write_to, minuend_hdul, file_path=file_path, comment_string=comment_string)

    return difference


# regex pattern for finding everything between brackets '[]'
# filename1 = '/home/lee/Documents/bias_frames/c4d_170331_202028_zri.fits.fz'
# filename2 = '/home/lee/Documents/bias_frames/c4d_170331_191113_zri.fits.fz'
source_directory = '/home/lee/Documents'
filename1 = 'k4m_160531_050920_ori.fits.fz'
filename2 = 'k4m_161228_132947_dri.fits.fz'
biasframe_file = 'testaverage.fits.fz'
# Open a new ds9 instance, or if already open, access it
display = pyds9.DS9(target='display', start='-title display')
ds9 = pyds9.DS9(target='ds9')
ds9.set('file ' + filename1)
with fits.open(filename1) as minuend_hdul, fits.open(filename2) as subtrahend_hdul:
    data_list = frame_subtract(minuend_hdul, subtrahend_hdul)

    difference = frame_subtract(ds9, filename2)
    # for data in data_list:
    #     display.set('frame new')
    #     display.set_np2arr(data)

frame_subtract(filename1, filename2, file_path=source_directory, display_in_ds9=True, write_to='test_result.fits.fz')

test_result = fits.open('/home/lee/Documents/test_result.fits.fz')
test_result
test_result.close()

__ = frame_subtract(filename1, filename2, file_path=source_directory,
                    display_in_ds9 = True, write_to = 'test_result.fits.fz')