from ccd_tools import *

path = '/home/lee/Data/Gain2'
filename_list = get_filenames(path=path)


def _sigma_clipped_frame_stats(image_stack, **kwargs):
    """Notes on this method:

    It does not account for the area of the guassian curve that is removed by
    the clipping. This means that the std dev result is potentially off by a
    factor equal to the clipped area.

    However, things like cosmic rays are *false signal*, and by def do not
    contribute to the *real* signal. Therefore, clipping them out does not affect
    the std dev: or more accurately, clipping them affects the std dev greatly,
    by removing the spurious signal and moving the std towards the std dev of the
    *true* signal."""
    # transpose the lists of data, then take the mean. np.squeeze removes the leftover axis
    # returns a list to account for multiple frames
    clipped_stats_frame_list = [sigma_clipped_stats(np.stack([hdu_data for hdu_data in data_tuple]), axis=0, **kwargs) for data_tuple in zip(*image_stack)]

    # remove the leftover axis from the data results.
    # also, unpack the results in separate lists, for consistency with the other stats functions
    clipped_mean_frame_list = [np.squeeze(frame_stats[0]) for frame_stats in clipped_stats_frame_list]
    clipped_median_frame_list = [np.squeeze(frame_stats[1]) for frame_stats in clipped_stats_frame_list]
    clipped_stddev_frame_list = [np.squeeze(frame_stats[2]) for frame_stats in clipped_stats_frame_list]

    return clipped_mean_frame_list, clipped_median_frame_list, clipped_stddev_frame_list


def sigma_clipped_frame_stats(filename_list, path='.', writeto_filename=None, overwrite=False, **kwargs):
    """

    Parameters
    ----------
    filename_list
    path
    writeto_filename
    overwrite
    kwargs

    Returns
    -------

    """
    # unpack the data, ignoring None values
    with ExitStack() as fits_stack:
        hdul_list = [fits_stack.enter_context(fits.open(path + '/' + fits_name)) for fits_name in filename_list]
        image_stack = [[hdu.data for hdu in hdul if hdu.data is not None] for hdul in hdul_list]

    if writeto_filename is not None:
        clipped_mean_frame_list, clipped_median_frame_list, clipped_stddev_frame_list = _sigma_clipped_frame_stats(image_stack, **kwargs)
        comment_string = 'Changed data to the average of ' + str(len(filename_list)) + \
                         ' zero frames, of which this is the first'

        _write_average_data_to_file(clipped_mean_frame_list, writeto_filename, source_filename_list=filename_list,
                                    file_path=path, overwrite=overwrite, comment_string=comment_string)
        return clipped_mean_frame_list, clipped_median_frame_list, clipped_stddev_frame_list
    else:
        return _sigma_clipped_frame_stats(image_stack, **kwargs)

