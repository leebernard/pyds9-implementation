"""
This script takes illumination data from the Andor iXon Ultra 888
camera, and measures the variance as a function of signal level.

The data was taken at -60 C, in a partially illuminated room, through
a Nikkor 50mm lens with several uncharacterized filters on top. The
filters both reduce illumination, and make the illumination more even.
"""

from ccd_tools import *
from astropy.stats import sigma_clip


# retrieve everything from the bias directory, ignoring files that are not fits
main_path = '/home/lee/Data/illumination_data_3'

bias_path = main_path + '/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

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
# open one of the files. the filename list is semi-random, so just grab the first entry
ds9.set('fits '+ sub_dir_filenames[0][0])
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
            selected_data = hdul[0].data[selection.xmin:selection.xmax, selection.ymin:selection.ymax]
            print(np.mean(selected_data))

variance = []
signal = []
gain = []
for exposure in sub_dir_filenames:
    with fits.open(exposure[0]) as frame1:
        with fits.open(exposure[1]) as frame2:
            # bias subtract the frame
            frame1_data = frame1[0].data.astype('float32') - master_bias_frame
            frame2_data = frame2[0].data.astype('float32') - master_bias_frame

            # crop the data down to what was selected, and sigma clip it
            frame1_data = sigma_clip(frame1_data[selection.xmin:selection.xmax, selection.ymin:selection.ymax], sigma=5.0)
            frame2_data = sigma_clip(frame2_data[selection.xmin:selection.xmax, selection.ymin:selection.ymax], sigma=5.0)

            frame_diff = frame1_data- frame2_data  # order of subtraction is arbitrary
            frame_var = np.var(frame_diff)
            frame_signal = np.median(np.asarray([frame1_data, frame2_data]))

            # diagnostic stuff
            print('differential average:', np.mean(frame_diff))
            # print('fraction of signal:', np.mean(frame_diff)/frame_signal)
            print('expected error:', np.sqrt(frame_var/frame_diff.size))
            # # display the frame difference
            # display_data(frame_diff)
            # histogram the frame difference
            # plt.figure(exposure[0])
            # plt.hist(frame_diff.flatten(), bins=90)
            # plt.figure(exposure[0] + 'signal')
            # plt.hist(frame1_data.flatten(), bins=100)
            # plt.show()

            # store signal
            signal.append(frame_signal)
            # store variance
            variance.append(frame_var)
            # store calculated gain, gain=signal/var * sqrt(2)
            gain.append(frame_signal/frame_var * np.sqrt(2))

plt.figure('gain vs signal')
plt.scatter(signal, gain)

plt.figure('ptc')
plt.scatter(signal, variance)



