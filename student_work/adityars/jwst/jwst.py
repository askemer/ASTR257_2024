import pandas as pd
from astropy import units as u
from astropy import constants as c

def pipe(inp, *funcs):
    r = inp
    for f in funcs:
        r = f(r)

    return r

def load_model_spectrum(fname):
     s = pd.read_csv(fname, header=1, sep="\s\s+", engine="python")
     s = s.rename(columns={"microns" : "Wavelength (um)", "Flux (erg/cm^2/s/Hz)" : "Flux (Jy)"})
     return s

def rescale_to_jwst_flux(spectrum):
    scaling = ((c.R_jup / (5.7 * u.pc)) ** 2).decompose() * 1000
    spectrum["Flux (Jy)"] *= scaling
    spectrum = spectrum.sort_values("Wavelength (um)")
    spectrum = spectrum.rename(columns={"Flux (Jy)" : "Flux (mJy)"})
    return spectrum

def drop_points(spectrum):
    return spectrum.iloc[::2]

def save_for_etc(spectrum, fname="sp_jwst_etc.dat"):
    spectrum.to_csv(fname, index=False, sep=" ", header=False)

if __name__ == "__main__":
    pipe("sp_t350g100nc_m0.0", load_model_spectrum, rescale_to_jwst_flux, drop_points, save_for_etc)
