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


folder = "/home/ev/UCSC/Classes/ASTR_257/ASTR257_2022/Project 1--Pluto/DiazDavisData/"
for i,file in enumerate(sorted(os.listdir(folder))):
    if not file.startswith('.'):
        data = fits.getdata(folder+file)
        header = fits.getheader(folder+file)

        fig = plt.figure(figsize=(3,3))
        print(file)
        # plot the raw image
        plt.title(file)
        plt.imshow(data, cmap = 'Blues_r', vmin=900, vmax = 1050)
        plt.colorbar();
