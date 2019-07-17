"""
A script for performing dark current analysis on the Andor iXon Ultra 888 camera.

This will bias-subtract and calculate dark current, at various temperatures.
Of particular note is how the CCD exhibits two dark currents: a low dark current
in most pixels, and a higher dark current in about 1% of pixels. Both dark currents
should be characterized as a function of temperature.
"""

from ccd_tools import *


'''
Starting out with -25 C
'''

# retrieve everything from the bias directory, ignoring files that are not fits
bias_path = '/home/lee/Data/darkcurrent_25c/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

bias_mean, bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(bias_files, path=bias_path, sigma=4.0)

bias_25c = bias_mean[0]

data_path = '/home/lee/Data/darkcurrent_25c/'

