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
how about -40c
'''

# use the gain published by Andor for 1MHz, 3.32 e-/ADC
gain = 3.32

# retrieve everything from the bias directory, ignoring files that are not fits
bias_path = '/home/lee/Data/darkcurrent_40c/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

# calculate a master bias for this temperature
bias_mean, bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(bias_files, path=bias_path, sigma=4.0)
bias_40c = bias_mean[0]


# pull the 480 data
path_480 = '/home/lee/Data/darkcurrent_40c/480s_exposure'
filenames_480 = get_filenames(path_480, extension='.fits')

# take the average of the frames, in counts
darkcurrent_480_average, _, _ = sigma_clipped_frame_stats(filenames_480, path=path_480, sigma=4.0)

# bias subtract, and then transform to e- from adc
darkcurrent_480_average = (darkcurrent_480_average - bias_40c)*gain

# show the histogram of the average dark current
plt.figure('dark current histrogram at 480s, -40c')
plt.hist(darkcurrent_480_average.flatten(), bins=120, range=(darkcurrent_480_average.min(), 2000))

# now sigma clip the data, to show the primary
primary_darkcurrent = sigma_clip(darkcurrent_480_average, sigma=4.0)
plt.figure('primary dark current, 480s, -40c')
plt.hist(primary_darkcurrent.compressed(), bins=120, range=(darkcurrent_480_average.min(), 2000))

# flip the mask, and histogram it
leftovers = np.ma.array(primary_darkcurrent.data, mask=~primary_darkcurrent.mask)
plt.figure('what is left over after removing the primary dark current pixels')
plt.hist(leftovers.compressed(),  bins=120, range=(darkcurrent_480_average.min(), 2000))

# clip out the extra stuff.
secondary_darkcurrent = sigma_clip(leftovers, sigma=2.0)
plt.figure('secondary dark current, 480s, -40c')
plt.hist(secondary_darkcurrent.compressed(),  bins=120, range=(darkcurrent_480_average.min(), 2000))

plt.show()

print('average taken at -40c, 480s exposure')
print('primary dark current:', np.mean(primary_darkcurrent)/480)
print('secondary dark current:', np.mean(secondary_darkcurrent)/480)

# ratio of primary dc pixels to secondary dc pixels
# assume area under the gaussian curve is equal to number of pixels
# correct for the area that has been clipped out in the secondary dark current
# 2 sigma --> .954 of total area
secondary_pixels_total = secondary_darkcurrent.count()/.954
# for the primary dc, sigma is 4 --> about 1.00 of total area
primary_pixels_total = primary_darkcurrent.count()

print('fraction of pixels that have primary dc value:', primary_pixels_total/secondary_darkcurrent.size)
print('fraction of pixels that has secondary dc value:', secondary_pixels_total/secondary_darkcurrent.size)

