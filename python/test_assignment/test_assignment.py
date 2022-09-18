# coding: utf-8

'''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
|                 C h r i s t o p h e r   E v a n   D a v i s                 |
|                                                                             |
|                               A S T R   2 5 7                               |
|                                                                             |
|    D e p t .   o f   A s t r o n o m y   a n d   A s t r o p h y s i c s    |
|                                                                             |
|    U n i v e r s i t y   o f   C a l i f o r n i a   S a n t a   C r u z    |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
'''

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["savefig.dpi"] = 300
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
import os, sys, shutil, subprocess
# import smart
import multiprocessing
import platform
import time
import colorcet as cc
import math
# import photochem_utils as pu
from decimal import Decimal
# import pubchempy as pcp
rcParams.update({'font.family':'sans-serif'})
import astropy as ap
import astropy.io.fits as fits
import photutils


# In[ ]:


# load image
img = fits.getdata('/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/python/test_assignment/test.fits', ext = 0)
print(img.shape)


# In[ ]:


# plot the raw image
plt.imshow(img, cmap = 'Blues_r');


# In[ ]:


# back up the original image
bkp_img = img
# set NaNs to zero
nan_mask = [[np.isnan(img[i,j]) for j in range(len(img[0]))] for i in range(len(img[0]))]
img[nan_mask] = 0

# plot the image w/ NaNs = 0
plt.imshow(img, cmap = 'Blues_r');


# In[ ]:


# find the brightest pixel
x_bright, y_bright = np.where(img == np.max(img))


# save the 40x40 region around the brightest pixel as a new fits file
hdu = fits.PrimaryHDU(img[int(x_bright)-19:int(x_bright)+20, int(y_bright)-19:int(y_bright)+20])
hdul = fits.HDUList(hdu)
hdul.writeto('/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/python/test_assignment/test2.fits', overwrite = True)


# load in the fits file and plot to verify it worked
img2 = fits.getdata('/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/python/test_assignment/test2.fits', ext = 0)
plt.imshow(img2, cmap = 'Blues_r');


# In[ ]:


# measure the centroid of this star, using the brightest pixel as a guess,
# and plot it to verify
cntrd = photutils.centroids.centroid_sources(img, x_bright, y_bright, box_size = 41)
plt.imshow(img, cmap = 'Blues_r')
plt.plot(cntrd[0], cntrd[1], c = "xkcd:red", marker = '.')

print(cntrd)
