"""
This script takes illumination data from the Andor iXon Ultra 888
camera, and measures the variance as a function of signal level.

The data was taken at -60 C, in a partially illuminated room, through
several uncharacterized filters. This was taken without a lens. The
filters both reduce illumination, and make the illumination more even.
"""

from ccd_tools import *
from astropy.stats import sigma_clip
from scipy.stats import norm
from scipy.stats import linregress


def run_gain_analysis(main_path):
    # retrieve everything from the bias directory, ignoring files that are not fits
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

    # start of main calculations:
    sub_dir_list = get_filenames(main_path, extension='exposure', include_path=True)
    print(sub_dir_list)  # print to check the output

    # this makes a list of lists, with each entry in the outer list corresponding to a
    # subdirectory, and each inner list being the file names withing that subdirectory
    sub_dir_filenames = []
    for sub_dir in sub_dir_list:
        # retrieve the filenames, checking the extension to make sure
        sub_dir_filenames.append(get_filenames(sub_dir, extension='.fits', include_path=True))

    # open an instance of DS9, to select a fairly flat region of data
    pyds9.ds9_targets()
    ds9 = pyds9.DS9(target='display')

    # in DS9, make a box region that avoids the edges of the image, to avoid any sort of trail-off or other edge
    # weirdness. Make sure that this region is flat to within a factor of two
    input('Pause while you select data. Press enter to continue')
    selection = get_ds9_region(ds9, get_data=False)

    # now pull the signal analysis: variance, signal level, and calculated gain
    variance = []
    signal = []
    gain = []
    exposure_time = []
    for exposure in sub_dir_filenames:
        for n, file in enumerate(exposure):
            with fits.open(exposure[n]) as frame1:
                with fits.open(exposure[n - 1]) as frame2:
                    # store the exposure time
                    time = frame1[0].header['exposure']
                    exposure_time.append(float(time))

                    # bias subtract the frame
                    frame1_data = frame1[0].data.astype('float32') - master_bias_frame
                    frame2_data = frame2[0].data.astype('float32') - master_bias_frame

                    # crop the data down to what was selected, and sigma clip it
                    frame1_data = sigma_clip(frame1_data[selection.ymin:selection.ymax, selection.xmin:selection.xmax],
                                             sigma=5.0)
                    frame2_data = sigma_clip(frame2_data[selection.ymin:selection.ymax, selection.xmin:selection.xmax],
                                             sigma=5.0)

                    frame_diff = frame1_data - frame2_data  # order of subtraction is arbitrary
                    frame_var = np.var(frame_diff)
                    frame_signal = np.median(np.asarray([frame1_data, frame2_data]))

                    # diagnostic stuff
                    # print('differential average:', np.mean(frame_diff))
                    # print('fraction of signal:', np.mean(frame_diff)/frame_signal)
                    # print('expected error:', np.sqrt(frame_var/frame_diff.size))
                    # # display the frame difference
                    # display_data(frame_diff)

                    # # histogram the frame difference
                    # plt.figure(exposure[n])
                    # plt.hist(frame_diff.compressed(), bins=np.arange(frame_diff.min(), frame_diff.max(), step=1), density=True)
                    # # add the fit of a guassian to this
                    # gmean, gstd = norm.fit(frame_diff.compressed())
                    # print('norm fit parameters:', gmean, gstd)
                    # print('std dev of data used in fit:', np.std(frame_diff.compressed()))
                    # xmin, xmax = plt.xlim()
                    # x = np.linspace(xmin, xmax, 1000)
                    # y = norm.pdf(x, gmean, gstd)
                    # plt.plot(x, y)

                    # # histogram one of the frames
                    # plt.figure(exposure[n] + 'signal')
                    # plt.hist(frame1_data.flatten(), bins=np.arange(-350.5, 350.5, step=1) + 3800)
                    # plt.show()

                    # store signal
                    signal.append(frame_signal)
                    # store variance
                    variance.append(frame_var)
                    # store calculated gain, gain=signal/var * sqrt(2)
                    gain.append(frame_signal / frame_var * 2)

    plt.figure('gain vs signal')
    plt.scatter(signal, gain)
    plt.ylim(0, 4)
    plt.title('gain measured from difference frames, as a function of signal from those frames.')
    plt.ylabel('gain (e-/ADC)')
    plt.xlabel('signal from a difference of frames (ADC)')

    plt.figure('ptc')
    plt.scatter(signal, variance)
    plt.title('photon transfer curve')
    plt.xlabel('Signal (ADC)')
    plt.ylabel('Variance (ADC)')

    # clipped the signal for where the gain deviates. Gain should be about the same for all values
    clipped_gain = sigma_clip(gain, sigma=3.0)
    clipped_signal = np.ma.array(signal, mask=clipped_gain.mask)
    fitresult = linregress(clipped_signal.compressed(), np.ma.array(variance, mask=clipped_signal.mask).compressed())
    slope = fitresult[0]
    intercept = fitresult[1]
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 1000)
    y = intercept + slope * x
    plt.plot(x, y, label='fit to slope')
    plt.legend()

    measured_gain = 1 / slope * 2
    print('Gain measured from slope of ptc:', measured_gain)

    plt.figure('exposure time vs signal')
    plt.scatter(exposure_time, np.asarray(signal) * measured_gain)
    plt.title('signal as a function of integration time')
    plt.xlabel('itegration time (s)')
    plt.xscale('log', basex=2)
    plt.ylabel('signal (e-)')
    plt.yscale('log')

    return measured_gain


# retrieve everything from the bias directory, ignoring files that are not fits
data_path = '/home/lee/Data/illumination_data_gain1_1MHz'
# data_path = '/home/lee/Data/illumination_data_gain2_100kHz'

gain = run_gain_analysis(data_path)


# # check the average of the selected region
# # while I'm at it, extract the data into a list
# region_data = []
# for filename_list in sub_dir_filenames:
#     print('-----------------------------------------------------')
#     print(filename_list)
#     print('Naive Average')
#     for file in filename_list:
#         with fits.open(file) as hdul:
#             # print(np.mean(hdul[0].data))
#             data_frame = hdul[0].data - master_bias_frame
#             selected_data = data_frame[selection.ymin:selection.ymax, selection.xmin:selection.xmax]
#             region_data.append(selected_data)
#             print(np.mean(selected_data))
#
# # check the data for patterns
# region_stack = np.stack(region_data)
# median_region = np.median(region_stack, axis=0)
# plt.figure('median of the region data')
# plt.hist(median_region.flatten(), bins=np.arange(-300.5, 300.5, step=1) + 3801)
# plt.title('Median of the region stack\n binned from 3500.5 to 4100.5')
# plt.figure('median of the region data, different bins')
# plt.hist(median_region.flatten(), bins=np.arange(-300, 300, step=1) + 3800)
# plt.title('Median of the region stack\n binned from 3500 to 4100')
# plt.figure('single region data')
# plt.hist(region_data[0].flatten(), bins=np.arange(-300, 300, step=1) + 3800)
# plt.title('One of the regions, pulled from the stack\n binned from 3500 to 4100')
