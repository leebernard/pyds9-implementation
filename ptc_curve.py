"""
This script takes illumination data from the Andor iXon Ultra 888
camera, and measures the variance as a function of signal level.

The data was taken at -60 C, in a partially illuminated room, through
a Nikkor 50mm lens with several uncharacterized filters on top. The
filters both reduce illumination, and make the illumination more even.
"""

from ccd_tools import *
from astropy.stats import sigma_clip

'''
how about -40c
'''

# retrieve everything from the bias directory, ignoring files that are not fits
main_path = '/home/lee/Data/illumination_data'

bias_path = main_path + '/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

# calculate a master bias for this temperature
bias_mean, bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(bias_files, path=bias_path, sigma=4.0)
master_bias_frame = bias_mean[0]



