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
import astropy.visualization as apv
import astropy.io.fits as fits
import photutils


# In[ ]:


'''
reduce the data for the first observing night
'''

# make empty numpy arrays of shape (5,) for our different frames
biases   = np.zeros(5, dtype = object)
darks    = np.zeros(5, dtype = object)
flats    = np.zeros(5, dtype = object)
sciences0 = np.zeros(5, dtype = object)

folder = "/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto0/"

for file in sorted(os.listdir(folder)):
    if file.endswith('.fits') and 'RSF' not in file:
        i = int(file.split('0')[1].replace('.fits',''))-1
        data = fits.getdata(folder+file)
        header = fits.getheader(folder+file)

        print(file,header['EXPTIME'])

        # add the files to respective arrays
        if 'bias' in file:
            biases[i] = data
        elif 'dark' in file:
            darks[i] = data
        elif 'flat' in file:
            flats[i] = data
        elif 'science' in file:
            sciences0[i] = data

    # # plot the data for a sanity check
    # zsi = apv.ZScaleInterval()
    # v = zsi.get_limits(data)
    # fig = plt.figure(figsize=(5,5))
    # plt.title(file)
    # plt.imshow(data, cmap = 'Blues_r', vmin=v[0], vmax = v[1])
    # plt.colorbar()
    # plt.tight_layout()
    # plt.savefig("/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto0/{:s}".format(file.replace('.fits','.png')))
    # plt.close();

# stack arrays so that the are 3d
biases    = np.stack(biases)
darks     = np.stack(darks)
flats     = np.stack(flats)
sciences0 = np.stack(sciences0)

# convert all the values to floats so arithmetic works
biases    = biases.astype('float64')
darks     = darks.astype('float64')
flats     = flats.astype('float64')
sciences0 = sciences0.astype('float64')

# take the mean of each type of frame
m_biases    = np.median(biases, axis = 0)
m_darks     = np.median(darks, axis = 0)

# normalize each flat
for i,flat in enumerate(flats):
    flats[i]  = flat - m_darks # subtracting median dark (has bias in it)
    norm_flat = np.median(flats[i]) # find median of each flat
    flats[i]  = flats[i]/norm_flat # divide each flat by its own meadian

m_flats     = np.median(flats, axis = 0)
m_sciences0 = np.median(sciences0, axis = 0)

# subtract the biases from the other frames
dark     = m_darks - m_biases
flat     = m_flats - m_biases
science0 = m_sciences0 - m_biases

# subtract the dark field from the flat field, accounting for exp. times
flat = flat - dark/3.

# normalize the flat field
norm_flat = np.median(flat)
flat      = flat/norm_flat

# reduce science frame using the calibration frames
rsf0 = (science0 - dark)/flat

# plot
zsi = apv.ZScaleInterval()
v = zsi.get_limits(rsf0)
fig = plt.figure(figsize=(5,5))
plt.title('Reduced Science Frame --- Pluto 0')
plt.imshow(rsf0, cmap = 'Blues_r', vmin=v[0], vmax = v[1])
plt.colorbar()
plt.tight_layout()
plt.savefig("/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto0/RSF_Pluto0")
plt.close();

# save the reduced image as a fits file
hdu = fits.PrimaryHDU(rsf0)
hdul = fits.HDUList(hdu)
hdul.writeto('/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto0/RSF_Pluto0.fits', overwrite = True)


# In[ ]:


'''
repeat the reduction for the second observing night
'''

# make empty numpy arrays of shape (5,) for our science frames
sciences1 = np.zeros(5, dtype = object)

folder = "/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto1/"

for file in sorted(os.listdir(folder)):
    if file.endswith('.fits'):
        i = int(file.replace('science','').replace('.fits',''))-1
        data = fits.getdata(folder+file)
        header = fits.getheader(folder+file)

        print(file,header['EXPTIME'])

        sciences1[i] = data


# take the mean of each type of frame
m_sciences1 = np.median(sciences1, axis = 0)

# convert all the values to floats so arithmetic works
m_sciences1 = m_sciences1.astype('float64')

# subtract the biases from the other frames
science1 = m_sciences1 - m_biases

# create the reduced science frame from the other medianed frames
rsf1 = (science1 - dark)/flat

# plot
zsi = apv.ZScaleInterval()
v = zsi.get_limits(rsf1)
fig = plt.figure(figsize=(5,5))
plt.title('Reduced Science Frame --- Pluto 1')
plt.imshow(rsf1, cmap = 'Blues_r', vmin=v[0], vmax = v[1])
plt.colorbar()
plt.tight_layout()
plt.savefig("/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto1/RSF_Pluto1")
plt.close();

# save the reduced image as a fits file
hdu = fits.PrimaryHDU(rsf1)
hdul = fits.HDUList(hdu)
hdul.writeto('/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/ArcysData/pluto1/RSF_Pluto1.fits', overwrite = True)
