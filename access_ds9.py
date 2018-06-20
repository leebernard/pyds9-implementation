# test space for learning to use pyds9

import pyds9 as ds9

# show targets
print(ds9.ds9_targets())

# need to have an instance of ds9 running
# define the ds9 object.  Calls class DS9
d = ds9.DS9()

# show the open file in the target instance of ds9
d.get('file')


# show the xpa id of target instance
d.access()

# 