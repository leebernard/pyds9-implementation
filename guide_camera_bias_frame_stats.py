from ccd_tools import *

# 0.1MHz Gain 2: least sensitive gain, slowest read speed
path = '/home/lee/Data/Gain2'
filenames = get_filenames(path=path)

# do stats on the first 5 files
bias_stats = []
bias_clipped_stats = []
for file in filenames[:5]:
   with fits.open(path+'/'+file) as hdul:
      bias_stats.append(image_stats(hdul[0].data))
      bias_clipped_stats.append(image_stats(hdul[0].data, sigma_clip=True, sigma=5.0))

stacked_bias_mean, stacked_bias_median, stacked_bias_stddev = sigma_clipped_frame_stats(filenames, writeto_filename='test', path=path, sigma=4.0)

# Do the same for 1MHz, Gain 1:
