

import numpy as np

from scipy import integrate


class empty: pass



def bin_width(log10ages, log10Zs):

    """ The age width of each bin """

    BW = np.zeros((len(log10ages), len(log10Zs)))
    min_age = 0
    for ia, log10age in enumerate(log10ages[:-1]):
        max_age = int(10**np.mean([log10ages[ia+1],log10age])) # years
        BW[ia,:] = max_age-min_age
        min_age = max_age
    return BW


def single_metallicity(log10ages, log10Zs, func, p):

    SFZH = np.zeros((len(log10ages), len(log10Zs)))

    # --- determine the metallicity bin in which to place all the SF

    iZ = (np.abs(log10Zs - p['log10Z'])).argmin()

    min_age = 0
    for ia, log10age in enumerate(log10ages[:-1]):
        max_age = int(10**np.mean([log10ages[ia+1],log10age])) # years
        sf = integrate.quad(func, min_age, max_age, args=(p))[0]
        SFZH[ia,iZ] = sf
        min_age = max_age

    # --- normalise
    SFZH /= np.sum(SFZH)
    SFZH *= 10**p['log10M*']
    SFR = None

    return SFZH, SFR




# --- constant star formation

def constant_func(t, p):
    if t < 10**(p['log10_duration']):
        return 1
    else:
        return 0

def constant(log10ages, log10Zs, p):
    return single_metallicity(log10ages, log10Zs, constant_func, p)



# --- constant star formation

def exponential_func(t, p):
    p['tau'] = 1E6*p['tau_Myr']
    if t < 10**(p['log10_duration']):
        return np.exp(-t/p['tau'])
    else:
        return 0

def exponential(log10ages, log10Zs, p):
    return single_metallicity(log10ages, log10Zs, exponential_func, p)





# --- quenched star formation

def quenched_func(t, p):
    if t < 10**(p['log10_total_duration']):
        if t > 10**(p['log10_quenched_duration']):
            return 1
        else:
            return 0
    else:
        return 0

def quenched(log10ages, log10Zs, p):
    return single_metallicity(log10ages, log10Zs, quenched_func, p)













# def CSFH(SPS, p, redshift=0.0, log10Z=-2.):
#
#     log10ages = SPS.grid['log10age']
#
#     iZ = (np.abs(SPS.grid['log10Z'] - log10Z)).argmin()
#
#     SFZH = np.zeros((len(SPS.grid['log10age']), len(SPS.grid['log10Z'])))
#
#     from astropy.cosmology import WMAP9 as cosmo
#
#     dz = 0.01
#     z = np.arange(10., redshift, -dz)
#
#
#
#     ages = -(cosmo.age(z).to('yr').value - cosmo.age(redshift).to('yr').value)
#
#     csfh = lambda z, p1, p2, p3, p4: p1 * (1+z)**p2 / (1+((1+z)/p3)**p4 )
#
#     sfrd = csfh(z, *p)
#
#     f = lambda x: np.interp(x, ages, sfrd[::-1])
#
#
#     start = 0.
#     for ia, log10age in enumerate(log10ages[:-1]):
#
#         end = 10**log10age
#
#         # determine the amount of star formation in this bin
#
#         sf = integrate.quad(f, start, end)[0]
#
#         SFZH[ia, iZ] = sf
#         # print('{0:.1f}-{1:.1f}: SF={2:.2f}'.format(start/1E6, end/1E6, sf))
#
#         start = end
#
#     SFZH /= np.sum(SFZH)
#
#     SFR = csfh(redshift, *p)
#
#     return SFZH, SFR







    # --- now rescale the mass




def instantaneous(log10ages, log10Zs, p):

    # --- returns the SFZH for an instantaneous SFH and single metallicity. At the moment this just chooses a single metallicity.

    SFZH = np.zeros((len(log10ages), len(log10Zs)))

    # --- determine the metallicity bin in which to place all the SF

    iZ = (np.abs(log10Zs - p['log10Z'])).argmin()

    ia = (np.abs(log10ages - p['log10age'])).argmin()

    SFZH[ia, iZ] = 10**p['log10M*']

    SFR = 0.0

    return SFZH, SFR
