# pyds9 Implementaion

The module pyds9 is a python interface for using SAOImage DS9. This is a wrapper and extension of 
pyds9, for use in the UCO Lick CCD lab. The main module is ccd_tools. The primary purpose of this 
module is for easy extraction of segments of image data from DS9 into a python environment, by 
detecting regions that have been selected by the user in DS9.  It also provides basic set of tools 
for statistical analysis. This is all very much a work in progress.

Some basic tools for statistical analysis have been added. There are now functions for finding 
average, median, and sigma clipped average frames. The frame average is found by stacking the data
frames, with each image forming a layer of the stack. If there are extensions, each header data unit 
gets it's own stack. It is presumed that header data units correspond to amplifiers, allowing for 
convenient stacking. The average is taken along the stack axis for each pixel. 

The median is found in a similar fashion to the average, and is intended to do the same thing, but
with some robustness against outlying pixel values, while having a shorter runtime than sigma 
clipping. 

The sigma clipped average is the most robust, and operates a little differently. It leverages 
Astropy's SigmaClip to remove any pixels with outlying values. This works by find the standard 
deviation (sigma) and median of the data on a image layer. Any pixels beyond a multiple of the sigma 
(say 5 sigma) on that particular layer are masked. This masking procedure can be iterated through a 
specified number of times, and the sigma to clip to can also be specified. Any pixels that are not 
masked are included in the average along the stack axis. See Astropy.stats Sigma_Clip for more 
advanced options. 

Because clipping removes pixels from the stack axis, it is possible for pixel average results to 
have been produced from fewer data points than the total number of frames being averaged. This has 
a couple of consequences. One, as pixels are clipped out, the error for that pixel increases (this 
is also why clipping is done on the image layer, instead of along the stack axis). Another is that 
some pixels can be clipped out entirely. If a pixel is clipped out entirely, it's data value is 
returned as zero. This can occur for hot pixels, and along the data edge. Both of those cases are 
handled by keeping track of how many data points are included for a particular pixel, in a data 
structure the same shape as the average result. Each entry in an array of this pixel count tracker
corresponds to a pixel in the image data, and contains a value equivalent to:

    (# of frames being averaged) - (# of times this particular pixel has been clipped) 
    = (number of data points used to compute average)

A robust way of subtracting frames has also been added, and is complex enough to need a little 
explanation. It is intended for subtracting bias and dark frames from image data, and displaying
the result in DS9. It can do this for frames on file, frames open in DS9, and a combination thereof.
This adds a wrinkle: in the case where header data unit (HDU) extensions are being used, DS9 only 
passes one extension at a time, the one currently selected by the user. If one frame is open in DS9, 
and the other on disk, the subtraction function will automatically acquire the appropriate HDU 
extension. If both frames being subtracted are open in DS9, they are simply subtracted. If both 
frames are on disk, the subtraction routine will iterate through all HDU extensions. 

The result of the frame sub

See http://hea-www.harvard.edu/RD/pyds9/ for more information on pyds9, including a link to source 
files on github.

Getting Started
---------------
#### Prerequisites

Python 3.6  
astropy  
numpy  
matplotlib  
scipy  
cython  
scikit-image

Optional:
photutils


#### Installing
Create an env with the requisite packages. Example:

    name@labcomp:~$ conda create --name env_name python=3.6
    name@labcomp:~$ conda activate env_name
    (env_name) name@labcomp:~$ conda install astropy numpy matplotlib scipy cython scikit-image
    (env_name) name@labcomp:~$ conda install -c astropy photutils
    (env_name) name@labcomp:~$ pip install pyds9

ccd_tools is available on github, at https://github.com/leebernard/pyds9-implementation/tree/distribute

One way of installing ccd_tools is to download directly:

go to  
https://raw.githubusercontent.com/leebernard/pyds9-implementation/distribute/ccd_tools.py  
Click (File or right-click) > Save Page As...  
Save to a folder named ccd_tools, in the directory of your choice.

Once the file has been downloaded, add this line to the .bash_profile (or .profile) file:

    export PYTHONPATH=/your/chosen/directory/ccd_tools:$PYTHONPATH

Alternately, you can enter the following code every time you run Python

    (env_name) name@labcomp:~$ python
    ... (startup messages) ...
    >>> import sys
    >>> sys.path.insert(0, '/your/chosen/directory/ccd_tools')
    >>> from ccd_tools import *

#### Examples
Retrieve a region from DS9, and show the statistics. This retrieves the bias subtracted data as 
well, if a bias region is available.

    >>> from ccd_tools import *
    # automatcially detect ds9 and retrieve a selected region
    >>> exampleregion = get_ds9_region()
    # Region file format: DS9 version 4.1
    global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    image
    Region definition:  ['box(1104.5,968.5,135,205,0)']
    Bias from  /home/lee/Documents/k4m_160531_050920_ori.fits.fz[im1]
    Bias Section is [2066:2112,18:2065]
    ----------------
    Bias statistics:
    Mean: 3357.89
    Median: 3358.00
    Std: 4.50
    Bias statistics after bias subtraction: 
    Mean: 0.00
    Median: 0.11
    Std: 4.50
    >>> exampleregion.stats()
    -----------------------
    Region Data Statistics:
    Min pixel value: 5367.00
    Max pixel value: 39805.00
    Mean: 5709.54
    Median: 5523.00
    Std: 1930.76
    ---------------------
    Bias Subtracted Data Statistics:
    Min pixel value: 2009.11
    Max pixel value: 36447.11
    Mean: 2351.64
    Median: 2165.11
    Std: 1930.76

To see the documentation of a particular function, enter print(foo.\_\_doc\_\_)

           
    >>> print(exampleregion.sky_subtract.__doc__)
    
            This method calculates the background subtracted data, and stores it in
            the sky_sub attribute.
    
            This is a wrapper for sky_subtract. If bias subtracted data is available, 
            it uses that. Otherwise, it uses data in the data attribute. If neither are
             available, raises an exception.
    
            Parameters
            ----------
            mask: numpy.ndarray (bool), optional
                A boolean mask with the same shape as data, where a True value
                indicates the corresponding element of data is masked. Masked pixels
                are excluded when computing the statistics.
            kwargs: dict, optional
                Optional keywords for sky_subtract
    
            Raises
            ------
            ValueError
                If both the data attribute and the bias_sub attribute are set to none.
    
            See Also
            --------
            sky_subtract: returns background subtracted data
            
    >>> 

Open a new instance of DS9. 
    
    # make a new, seperate instance of ds9
    >>> ds9display = pyds9.DS9(target='display', start='-title display')
    >>>

To show the currently available instances of DS9, enter
    
    >>> pyds9.ds9_targets()
    ['DS9:display 7f000001:37337', 'DS9:ds9 7f000001:44345']
To select a specific instance (say if you opened multiple instances before starting python), enter

    >>> ds9 = pyds9.DS9(target='display')
    >>> 
Alternately, you can specify the XPA address. This is functionally the same as the above:

    >>> ds9 = pyds9.DS9(target='7f000001:37337')
    >>> 
    
#### Notes on Using the Region Class
If you are coming from a procedural language like C, you are probably not familiar with classes. If 
you know a language like C++, you may be familiar with classes and objects, but there are some 
idiosyncrasies with how classes are implemented in Python that are good to know about. The easiest 
way to understand objects is to work through a use example. Therefore, this section will give a 
brief run-down of how the region class operates. 

In the above example, the data from SAOImage DS9 is stored in an instance of the region class. A 
class is a structure for convenient handling of related data and functions. In this case, the region 
class is a wrapper for the various types of data associated with a SAOImage DS9 region, such as the 
source file name, the location of the region in the frame, and of course the pixel values. More 
information about the region class can be found in the documentation.

Instances of a class, referred to as objects, are unique in that changes to one instance will not 
affect any other instances. In this case the object was defined as 'exampleregion'. 

    >>> from ccd_tools import *
    >>> exampleregion = get_ds9_region()
    /home/lee/PycharmProjects/pyds9-implementation/ccd_tools.py:592: UserWarning: "Keyword 'BIASSEC' not found." Unable to perform bias subtraction.
      warnings.warn(message, category=UserWarning)
    >>> 

The pixel value data from the region selected in SAOImage DS9 is stored in exampleregion, under the 
.data attribute. This data can be accessed by calling the data attribute:

    >>> exampleregion.data
    array([[2663, 2657, 2658, ..., 2631, 2661, 2660],
           [2643, 2655, 2635, ..., 2625, 2637, 2637],
           [2619, 2648, 2622, ..., 2631, 2642, 2630],
           ...,
           [2653, 2643, 2633, ..., 2651, 2625, 2626],
           [2633, 2661, 2626, ..., 2631, 2660, 2614],
           [2628, 2635, 2665, ..., 2631, 2661, 2633]], dtype=uint16)

The get_ds9_region() fuction will also attempt to provide bias-subtracted data. In this case, 
get_ds9_region() was unable to find the bias section of the frame containing the region, because 
the default keyword was not found in the header. If this data had been found, it would of been 
stored under .bias_sub.

    >>> exampleregion.bias_sub
    >>> print(exampleregion.bias_sub)
    None

As can be seen, built in attributes have a default value of None. 

Python has a very flexible implementation of objects. Unlike in C++, all attributes can be modified 
at any time by the user. This also means that the user can add new attributes at will. 

    >>> exampleregion.data = 'oops'
    >>> print(exampleregion.data)
    oops
    >>> # creating a new data attribute, containing a string
    ... exampleregion.note_to_self = 'The data attribute was accidentally overwritten'
    


#### Issues
Sometimes, if multiple instances of DS9 are closed, the XPA name server seems to close. If it shows 
an error similar to this
    
    ValueError: no active ds9 running for target: display
Try running 
    
        >>> pyds9.ds9_xpans()

See the documentation for pyds9 for more information, at http://hea-www.harvard.edu/RD/pyds9/

Lee Bernard