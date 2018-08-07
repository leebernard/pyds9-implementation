# pyds9 Implementaion

The module pyds9 is a python interface for using SAOImage DS9. This is a wrapper and extension of 
pyds9, for use in the UCO Lick CCD lab. The main module is ccd_tools. The primary purpose of this 
module is for easy extraction of segments of image data from DS9 into a python environment, by 
detecting regions that have been selected by the user in DS9.  It also provides basic set of tools 
for statistical analysis. These tools will be expanded on in the future.

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
    
#### Issues
Sometimes, if multiple instances of DS9 are closed, the XPA name server seems to close. If it shows 
an error similar to this
    
    ValueError: no active ds9 running for target: display
Try running 
    
        >>> pyds9.ds9_xpans()

See the documentation for pyds9 for more information, at http://hea-www.harvard.edu/RD/pyds9/

Lee Bernard