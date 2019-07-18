"""
A script for performing dark current analysis on the Andor iXon Ultra 888 camera.

This will bias-subtract and calculate dark current, at various temperatures.
Of particular note is how the CCD exhibits two dark currents: a low dark current
in most pixels, and a higher dark current in about 1% of pixels. Both dark currents
should be characterized as a function of temperature.
"""

from ccd_tools import *
from astropy.stats import sigma_clip

'''
Starting out with -25 C
'''

# retrieve everything from the bias directory, ignoring files that are not fits
bias_path = '/home/lee/Data/darkcurrent_25c/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

# calculate a master bias for this temperature
bias_mean, bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(bias_files, path=bias_path, sigma=4.0)
bias_25c = bias_mean[0]

# retrieve the dark current frame data, and relevant metadata
data_path = '/home/lee/Data/darkcurrent_25c'
darkcurrent_filenames = get_filenames(data_path, extension='.fits')

darkcurrent_frames = []
temperature = []
exposure_time = []
for filename in darkcurrent_filenames:
    with fits.open(data_path + '/' + filename) as hdulist:
        print('data type:', hdulist[0].data.dtype)
        darkcurrent_frames.append(hdulist[0].data)

        print('temperature:', hdulist[0].header['temp'])
        temperature.append(hdulist[0].header['temp'])

        print('exposure time:', hdulist[0].header['exposure'])
        exposure_time.append(hdulist[0].header['exposure'])

# subtract off the bias
darkcurrent_frames = [frame - bias_25c for frame in darkcurrent_frames]
# convert to e- from ADC
# use the gain published by Andor for 1MHz, 3.32 e-/ADC
gain = 3.32
darkcurrent_frames = [frame*gain for frame in darkcurrent_frames]

# calculate naive dark current rate, and plot it
darkcurrent_rate = [np.mean(foo[0])/foo[1] for foo in zip(darkcurrent_frames, exposure_time)]
plt.scatter(exposure_time, darkcurrent_rate)

# calculate the sigma clipped dark current rate, and plot it
clipped_darkcurrent_frames = [sigma_clip(frame, sigma=4.0) for frame in darkcurrent_frames]
clipped_darkcurrent_rate = [np.mean(foo[0])/foo[1] for foo in zip(clipped_darkcurrent_frames, exposure_time)]
plt.scatter(exposure_time, clipped_darkcurrent_rate)

plt.figure()
plt.hist(darkcurrent_frames[4].flatten(), bins=100, range=(darkcurrent_frames[4].min(), 2000))

plt.figure()
plt.hist(clipped_darkcurrent_frames[4].compressed(), bins=100, range=(darkcurrent_frames[4].min(), 2000))

