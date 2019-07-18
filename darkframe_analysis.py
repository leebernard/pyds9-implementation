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

# use the gain published by Andor for 1MHz, 3.32 e-/ADC
gain = 3.32

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
darkcurrent_frames = [frame*gain for frame in darkcurrent_frames]

# calculate naive dark current rate, and plot it
darkcurrent_rate = [np.mean(foo[0])/foo[1] for foo in zip(darkcurrent_frames, exposure_time)]
plt.figure('exposure vs darkcurrent rate')
plt.scatter(exposure_time, darkcurrent_rate)

# calculate the sigma clipped dark current rate, and plot it
clipped_darkcurrent_frames = [sigma_clip(frame, sigma=4.0) for frame in darkcurrent_frames]
clipped_darkcurrent_rate = [np.mean(foo[0])/foo[1] for foo in zip(clipped_darkcurrent_frames, exposure_time)]
plt.scatter(exposure_time, clipped_darkcurrent_rate)


'''histogram a sample 240 exposure, see if I can isolate the two different dark currents'''
'''plt.figure()
plt.hist(darkcurrent_frames[4].flatten(), bins=100, range=(darkcurrent_frames[4].min(), 2000))

plt.figure()
plt.hist(clipped_darkcurrent_frames[4].compressed(), bins=100, range=(darkcurrent_frames[4].min(), 2000))

# flip the mask, and histogram it
secondary_darkcurrent = np.ma.array(clipped_darkcurrent_frames[4].data, mask=~clipped_darkcurrent_frames[4].mask)
plt.figure()
plt.hist(secondary_darkcurrent.compressed(), bins=100, range=(darkcurrent_frames[4].min(), 2000))

# clip out the extra stuff.
clipped_secondary_darkcurrent = sigma_clip(secondary_darkcurrent, sigma=2.0)
plt.figure('secondary dark current')
plt.hist(clipped_secondary_darkcurrent.compressed(), bins=100, range=(darkcurrent_frames[4].min(), 2000))

print('Sample taken at -25c, 240s exposure')
print('primary dark current:', np.mean(clipped_darkcurrent_frames[4])/exposure_time[4])
print('secondary dark current:', np.mean(clipped_secondary_darkcurrent)/exposure_time[4])

# ratio of primary dc pixels to secondary dc pixels
# assume area under the gaussian curve is equal to number of pixels
# correct for the area that has been clipped out in the secondary dark current
# 2 sigma --> .954 of total area
secondary_pixels_total = clipped_secondary_darkcurrent.count()/.954
# for the primary dc, sigma is 4 --> about 1.00 of total area
primary_pixels_total = clipped_darkcurrent_frames[4].count()

print('fraction of pixels that have primary dc value:', primary_pixels_total/clipped_secondary_darkcurrent.size)
print('fraction of pixels that has secondary dc value:', secondary_pixels_total/clipped_secondary_darkcurrent.size)
'''

'''okay, lets do it again, but with the 5 frames stacked and averaged, to reduce noise'''

# pull the 240 data
path_240 = '/home/lee/Data/darkcurrent_25c/240s_exposure'
filenames_240 = get_filenames(path_240, extension='.fits')

# take the average of the frames, in counts
darkcurrent_240_average, _, _ = sigma_clipped_frame_stats(filenames_240, path=path_240, sigma=4.0)

# bias subtract, and then transform to e- from adc
darkcurrent_240_average = (darkcurrent_240_average - bias_25c)*gain

# show the histogram of the average dark current
plt.figure('dark current histrogram at 240s, -25c')
plt.hist(darkcurrent_240_average.flatten(), bins=120, range=(darkcurrent_240_average.min(), 2000))

# now sigma clip the data, to show the primary
primary_darkcurrent = sigma_clip(darkcurrent_240_average, sigma=4.0)
plt.figure('primary dark current, 240s, -25c')
plt.hist(primary_darkcurrent.compressed(), bins=120, range=(darkcurrent_240_average.min(), 2000))

# flip the mask, and histogram it
leftovers = np.ma.array(primary_darkcurrent.data, mask=~primary_darkcurrent.mask)
plt.figure('what is left over after removing the primary dark current pixels')
plt.hist(leftovers.compressed(),  bins=120, range=(darkcurrent_240_average.min(), 2000))

# clip out the extra stuff.
secondary_darkcurrent = sigma_clip(leftovers, sigma=2.0)
plt.figure('secondary dark current, 240s, -25c')
plt.hist(secondary_darkcurrent.compressed(),  bins=120, range=(darkcurrent_240_average.min(), 2000))

print('average taken at -25c, 240s exposure')
print('primary dark current:', np.mean(primary_darkcurrent)/240)
print('secondary dark current:', np.mean(secondary_darkcurrent)/240)

# ratio of primary dc pixels to secondary dc pixels
# assume area under the gaussian curve is equal to number of pixels
# correct for the area that has been clipped out in the secondary dark current
# 2 sigma --> .954 of total area
secondary_pixels_total = secondary_darkcurrent.count()/.954
# for the primary dc, sigma is 4 --> about 1.00 of total area
primary_pixels_total = primary_darkcurrent.count()

print('fraction of pixels that have primary dc value:', primary_pixels_total/secondary_darkcurrent.size)
print('fraction of pixels that has secondary dc value:', secondary_pixels_total/secondary_darkcurrent.size)


'''again, but with the 60s'''

path_60 = '/home/lee/Data/darkcurrent_25c/60s_exposure'
filenames_60 = get_filenames(path_60, extension='.fits')

# take the average of the frames, in counts
darkcurrent_60_average, _, _ = sigma_clipped_frame_stats(filenames_60, path=path_60, sigma=4.0)

# bias subtract, and then transform to e- from adc
darkcurrent_60_average = (darkcurrent_60_average - bias_25c)*gain

# show the histogram of the average dark current
plt.figure('dark current histrogram at 60s, -25c')
plt.hist(darkcurrent_60_average.flatten(), bins=120, range=(darkcurrent_60_average.min(), 2000))

# now sigma clip the data, to show the primary
primary_darkcurrent = sigma_clip(darkcurrent_60_average, sigma=4.0)
plt.figure('primary dark current, 60s, -25c')
plt.hist(primary_darkcurrent.compressed(), bins=120, range=(darkcurrent_60_average.min(), 2000))

# flip the mask, and histogram it
leftovers = np.ma.array(primary_darkcurrent.data, mask=~primary_darkcurrent.mask)
plt.figure('what is left over after removing the primary dark current pixels')
plt.hist(leftovers.compressed(),  bins=120, range=(darkcurrent_60_average.min(), 2000))

# clip out the extra stuff.
secondary_darkcurrent = sigma_clip(leftovers, sigma=2.0)
plt.figure('secondary dark current, 60s, -25c')
plt.hist(secondary_darkcurrent.compressed(),  bins=120, range=(darkcurrent_60_average.min(), 2000))

print('average taken at -25c, 60s exposure')
print('primary dark current:', np.mean(primary_darkcurrent)/60)
print('secondary dark current:', np.mean(secondary_darkcurrent)/60)

# ratio of primary dc pixels to secondary dc pixels
# assume area under the gaussian curve is equal to number of pixels
# correct for the area that has been clipped out in the secondary dark current
# 2 sigma --> .954 of total area
secondary_pixels_total = secondary_darkcurrent.count()/.954
# for the primary dc, sigma is 4 --> about 1.00 of total area
primary_pixels_total = primary_darkcurrent.count()

print('fraction of pixels that have primary dc value:', primary_pixels_total/secondary_darkcurrent.size)
print('fraction of pixels that has secondary dc value:', secondary_pixels_total/secondary_darkcurrent.size)

