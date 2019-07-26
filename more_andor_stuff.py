"""
This script takes illumination data from the Andor iXon Ultra 888
camera, and measures the variance as a function of signal level.

The data was taken at -60 C, in a partially illuminated room, through
several uncharacterized filters. This was taken without a lens. The
filters both reduce illumination, and make the illumination more even.
"""

from ccd_tools import *
from astropy.stats import sigma_clip


# retrieve everything from the bias directory, ignoring files that are not fits
main_path = '/home/lee/Data/illumination_data_4'

bias_path = main_path + '/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

# calculate read noise, from a couple of arbitrary frames
with fits.open(bias_path + '/' + bias_files[0]) as bias_hdul:
    bias_frame1 = bias_hdul[0].data
with fits.open(bias_path + '/' + bias_files[5]) as bias_hdul:
    bias_frame6 = bias_hdul[0].data

print('Stats on', bias_files[0])
print(sigma_clipped_stats(bias_frame1, sigma=4.0))
print('stats on', bias_files[5])
print(sigma_clipped_stats(bias_frame6, sigma=4.0))

# calculate a master bias for this temperature
bias_mean, bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(bias_files, path=bias_path, sigma=4.0)
master_bias_frame = bias_mean[0]

sub_dir_list = get_filenames(main_path, extension='exposure', include_path=True)
print(sub_dir_list)  # print to check the output