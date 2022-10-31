import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from os import path, listdir
from astropy.visualization import ZScaleInterval
from tqdm import tqdm, trange

datapath = "../../../Project 3--Neptune AO/"
getlims = ZScaleInterval().get_limits

def mediandiv_and_zerocorr(x):
    x_norm = x / np.nanmedian(x)
    # x_norm[np.where(x_norm == 0)] = 1
    return x_norm

def view_as_ds9(img, ax=None, **kwargs):
    vmin, vmax = getlims(img)
    if vmin > vmax:
        vmin, vmax = vmax, vmin
    if ax is None:
        plt.imshow(img, vmin=vmin, vmax=vmax, **kwargs)
    else:
        ax.imshow(img, vmin=vmin, vmax=vmax, **kwargs)

def process_flats():
    flats = {"H": [], "J": [], "K": []}
    flatnames = listdir(path.join(datapath, "skyflats"))
    for f in flatnames:
        flats[f[8]].append(fits.open(path.join(datapath, "skyflats", f))[0].data)


    medflats = {k : mediandiv_and_zerocorr(np.median([mediandiv_and_zerocorr(x) for x in flats[k]], axis=0)) for k in flats}
    for k in medflats:
        fits.writeto(path.join(datapath, "medflats", f"{k}.fits"), medflats[k], overwrite=True)

    return medflats

def retrieve_flats():
    medflats = {}
    for filtername in ["H", "J", "K"]:
        medflats[filtername] = fits.open(path.join(datapath, "medflats", f"{filtername}.fits"))[0].data
    return medflats

wherenan = lambda x: np.where(np.isnan(x))

def get_sciences():
    sciences = {"H": [], "J": [], "K": []}
    scinames = filter(lambda x: len(x) > 13, listdir(path.join(datapath, "science")))
    for f in scinames:
        sciences[f[8]].append(fits.open(path.join(datapath, "science", f))[0].data.astype(np.float64))

    return sciences


yl, yh, xl, xh = 250, 1250, 600, 1600
yc = (yl + yh) / 2 
xc = (xl + xh) / 2
r = xc - xl
xmesh, ymesh = np.meshgrid(np.arange(xl, xh), np.arange(yl, yh))
mask = (xmesh - xc) ** 2 + (ymesh - yc) ** 2 < r ** 2

def crop_to_ao(img):
    img = img[yl:yh, xl:xh]
    return img * mask

argmax = lambda x: np.where(x == np.max(x))

def renorm_to_find_neptune(img):
    img = n.crop_to_ao(img)
    cmin, cmax = n.getlims(img)
    return np.nan_to_num(
        np.maximum(
            cmin, np.minimum(
                cmax, img
            )
        )
    )[200:-200, 200:-200]