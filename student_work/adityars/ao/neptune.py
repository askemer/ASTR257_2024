import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from os import path, listdir
from astropy.visualization import ZScaleInterval
from tqdm import tqdm, trange

datapath = "../../../Project 3--Neptune AO/"
getlims = ZScaleInterval().get_limits

filternames = ["H", "J", "K"]

get_idxs = {
    "H" : [1, 2, 3, 4],
    "J" : [0, 1, 2, 3],
    "K" : [1, 2, 3, 4]
}

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
    
    for filtername in sciences:
        sciences[filtername] = [x for (i, x) in enumerate(sciences[filtername]) if i in get_idxs[filtername]] 
        # python how dare you make something so easy into an O(N^2) operation. I will never forgive you 

    return sciences

def divide_flats(sciences, flats):
    with np.errstate(divide='ignore', invalid='ignore'):
        for k in sciences:
            for i in range(len(sciences[k])):
                sciences[k][i] /= flats[k]
                sciences[k][i][wherenan(sciences[k][i])] = 0

    return sciences


yl, yh, xl, xh = 500, 1000, 950, 1450

argmax = lambda x: np.where(x == np.max(x))

def renorm(img):
    cmin, cmax = getlims(img)
    return np.nan_to_num(
        np.maximum(
            cmin, np.minimum(
                cmax, img
            )
        )
    )

def renorm_and_crop_all(sciences, do_renorm=True):
    if sciences["H"][0].shape[0] == 2048:
        for f in ["H", "J", "K"]:
            for i in range(4):
                img = sciences[f][i][yl:yh, xl:xh]
                if do_renorm:
                    img = renorm(img)
                sciences[f][i] = img
    
    return sciences

def plot_all_images(sciences, filtername):
    f, axs = plt.subplot_mosaic("01;23")
    for i in range(4):
        img = sciences[filtername][i]
        ax = axs[str(i)]
        ax.set_xticks([])
        ax.set_yticks([])
        ax.imshow(img, cmap="gray")
        #view_as_ds9(img, ax, cmap="gray")

    f.suptitle(f"Dithered images of Neptune, {filtername} filter")

