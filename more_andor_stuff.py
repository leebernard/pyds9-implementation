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

# this makes a list of lists, with each entry in the outer list corresponding to a
# subdirectory, and each inner list being the file names withing that subdirectory
sub_dir_filenames = []
for sub_dir in sub_dir_list:
    # retrieve the filenames, checking the extension to make sure
    sub_dir_filenames.append(get_filenames(sub_dir, extension='.fits', include_path=True))

# print to check the output
for filename_list in sub_dir_filenames:
    print(filename_list)


# open an instance of DS9, to select a fairly flat region of data
ds9 = pyds9.DS9(target='display')
pyds9.ds9_targets()

# in DS9, make a box region that avoids the edges of the image, to avoid any sort of trail-off or other edge
# weirdness. Make sure that this region is flat to within a factor of two
input('Pause while you select data. Press enter to continue')


selection = get_ds9_region(ds9, get_data=False)

# check the average of the selected region
for filename_list in sub_dir_filenames:
    print('-----------------------------------------------------')
    print(filename_list)
    print('Naive Average')
    for file in filename_list:
        with fits.open(file) as hdul:
            # print(np.mean(hdul[0].data))
            selected_data = hdul[0].data[selection.ymin:selection.ymax, selection.xmin:selection.xmax]
            print(np.mean(selected_data))


