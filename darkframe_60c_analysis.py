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
Here it is: -60 C
'''

# use the gain published by Andor for 1MHz, 0.790 e-/ADC
gain = 0.790

# retrieve everything from the bias directory, ignoring files that are not fits
bias_path = '/home/lee/Data/darkcurrent_60c/bias_frames'
bias_files = get_filenames(bias_path, extension='.fits')

for bias_file in bias_files:
    with fits.open(bias_path + '/' + bias_file) as bias_frame:
        print(bias_file)
        data = bias_frame[0].data
        print('min:', np.min(data))
        print('max:', np.max(data))
        print('std:', np.std(data))
# calculate a master bias for this temperature
bias_mean, bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(bias_files, path=bias_path, sigma=3.0)
bias_60c = bias_mean[0]


# pull the 5400 data
path_5400 = '/home/lee/Data/darkcurrent_60c/5400s_exposure'
filenames_5400 = get_filenames(path_5400, extension='.fits')

# print the averages of each frame
for filename in filenames_5400:
    with fits.open(path_5400 + '/' + filename) as hdul:
        print(filename)
        # hdul[0].data
        # print(np.median(hdul[0].data))
        print(np.mean(sigma_clip(hdul[0].data, sigma=5.0)))
# take the average of the frames, in counts
darkcurrent_5400_average, _, _ = sigma_clipped_frame_stats(filenames_5400, path=path_5400, sigma=4.0)

# bias subtract, and then transform to e- from adc
darkcurrent_5400_average = (darkcurrent_5400_average - bias_60c)*gain

# show the histogram of the average dark current
plt.figure('dark current histrogram at 5400s, -60c')
plt.hist(darkcurrent_5400_average.flatten(), bins=120, range=(darkcurrent_5400_average.min(), 1000))

# now sigma clip the data, to show the primary
primary_darkcurrent = sigma_clip(darkcurrent_5400_average, sigma=4.0)
plt.figure('primary dark current, 5400s, -60c')
plt.hist(primary_darkcurrent.compressed(), bins=120, range=(darkcurrent_5400_average.min(), 1000))

# flip the mask, and histogram it
leftovers = np.ma.array(primary_darkcurrent.data, mask=~primary_darkcurrent.mask)
plt.figure('what is left over after removing the primary dark current pixels')
plt.hist(leftovers.compressed(),  bins=120, range=(darkcurrent_5400_average.min(), 1000))

# clip out the extra stuff.
secondary_darkcurrent = sigma_clip(leftovers, sigma=2.0)
plt.figure('secondary dark current, 5400s, -60c')
plt.hist(secondary_darkcurrent.compressed(),  bins=120, range=(darkcurrent_5400_average.min(), 2000))



print('average taken at -60c, 5400s exposure')
print('primary dark current:', np.mean(primary_darkcurrent)/5400)
print('secondary dark current:', np.mean(secondary_darkcurrent)/5400)

# ratio of primary dc pixels to secondary dc pixels
# assume area under the gaussian curve is equal to number of pixels
# correct for the area that has been clipped out in the secondary dark current
# 2 sigma --> .954 of total area
secondary_pixels_total = secondary_darkcurrent.count()/.954
# for the primary dc, sigma is 4 --> about 1.00 of total area
primary_pixels_total = primary_darkcurrent.count()

print('fraction of pixels that have primary dc value:', primary_pixels_total/secondary_darkcurrent.size)
print('fraction of pixels that has secondary dc value:', secondary_pixels_total/secondary_darkcurrent.size)

plt.show()

