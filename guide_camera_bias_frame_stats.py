from ccd_tools import *

# 0.1MHz Gain 2: least sensitive gain, slowest read speed
path = '/home/lee/Data/Gain2'
# match start of string, then look ahead infinitely to see if 'master' is not in string
pattern = re.compile(r'^(?!.*(master))')
filenames = get_filenames(path=path, pattern=pattern, extension='fits')

# do stats on the first 5 files
bias_stats = []
bias_clipped_stats = []
for file in filenames[:5]:
    with fits.open(path+'/'+file) as hdul:
        print('data type:', hdul[0].data.dtype)
        bias_stats.append(image_stats(hdul[0].data))
        bias_clipped_stats.append(image_stats(hdul[0].data, sigma_clip=True, sigma=4.0))

# print the unclipped bias results
print('bias frame stats:')
print('mean    |median  |std dev   ')
for stattuple in bias_stats:
    print(f'{stattuple[0]:.3f}', f'|{stattuple[1]:.3f}', f'|{stattuple[2]:.3f}')

# print the bias results, this time the clipped stuff
print('sigma clipped bias frame stats:')
print('mean    |median  |std dev   ')
for stattuple in bias_clipped_stats:
    print(f'{stattuple[0]:.3f}', f'|{stattuple[1]:.3f}', f'|{stattuple[2]:.3f}')

stacked_bias_mean, stacked_bias_median, stacked_bias_stddev = \
    sigma_clipped_frame_stats(filenames, writeto_filename='bias_gain2_100kHz_20190628_master.fits', path=path, sigma=4.0)

print('mean pixel count of the master bias frame:', f'{np.mean(stacked_bias_mean):.3f}')
print('median pixel count of same:', f'{np.median(stacked_bias_mean):.3f}')
print('standard deviation of same:', f'{np.std(stacked_bias_mean):.3f}')
# Do the same for 1MHz, Gain 1:

# now it's time to measure the read noise
# start an instance of ds9, to display a bias frame
ds9 = pyds9.DS9(target='bias-display', start='-title display')
pyds9.ds9_targets()
# open one of the files. the filename list is semi-random, so just grab the first entry
ds9.set('fits '+path+'/'+filenames[0])
# in DS9, make a box region that avoids the edges of the image, to avoid any sort of trail-off or other edge
# weirdness. Also make sure no cosmic rays are in it. Then select the box.
input('Pause while you select data. Press enter to continue')

bias_data = get_ds9_region(ds9).data
# calculate statistics on the data, such that it prints the results
print('Statistics without clipping')
image_stats(bias_data, sigma_clip=False, verbose=True)
# readnoise is the std dev of the data, after sigma clipping
print('------------------------------------')
print('Statistics after clipping to 4-sigma')
_ , _, readnoise = image_stats(bias_data, sigma_clip=True, sigma=4.0, verbose=True)
