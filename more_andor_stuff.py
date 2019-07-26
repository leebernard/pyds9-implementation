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

# display a histogram of the master bias file
plt.figure('master bias frame')
plt.hist(master_bias_frame.flatten(), bins=100, range=(480, 600))

# show what is not part of the gaussian
clipped_master_bias = sigma_clip(master_bias_frame, sigma=4.0)
leftovers = np.ma.array(clipped_master_bias.data, mask=~clipped_master_bias.mask)
print('number of clipped pixels:', leftovers.compressed().size)

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


variance = []
signal = []
gain = []
for exposure in sub_dir_filenames:
    for n, file in enumerate(exposure):
        with fits.open(exposure[n]) as frame1:
            with fits.open(exposure[n-1]) as frame2:
                # bias subtract the frame
                frame1_data = frame1[0].data.astype('float32') - master_bias_frame
                frame2_data = frame2[0].data.astype('float32') - master_bias_frame

                # crop the data down to what was selected, and sigma clip it
                frame1_data = sigma_clip(frame1_data[selection.xmin:selection.xmax, selection.ymin:selection.ymax], sigma=5.0)
                frame2_data = sigma_clip(frame2_data[selection.xmin:selection.xmax, selection.ymin:selection.ymax], sigma=5.0)

                frame_diff = frame1_data - frame2_data  # order of subtraction is arbitrary
                frame_var = np.var(frame_diff)
                frame_signal = np.median(np.asarray([frame1_data, frame2_data]))

                # diagnostic stuff
                # print('differential average:', np.mean(frame_diff))
                # print('fraction of signal:', np.mean(frame_diff)/frame_signal)
                # print('expected error:', np.sqrt(frame_var/frame_diff.size))
                # # display the frame difference
                # display_data(frame_diff)
                # histogram the frame difference
                plt.figure(exposure[n])
                plt.hist(frame_diff.flatten(), bins=np.arange(-350.5, 350.5, step=1))
                plt.figure(exposure[n] + 'signal')
                plt.hist(frame1_data.flatten(), bins=np.arange(-350.5, 350.5, step=1) + 3800)
                plt.show()

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
