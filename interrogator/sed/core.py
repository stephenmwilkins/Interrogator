
import numpy as np
import scipy.integrate

from . import IGM

from flare.photom import *

h = 6.626E-34*1E7 # erg/Hz
c = 3.E8 #m/s


class sed():

    def __init__(self, lam, description = False):

        self.description = description

        self.lam = lam # \AA
        self.lnu = np.zeros(self.lam.shape) # luminosity ers/s/Hz
        self.nu = 3E8/(self.lam/1E10) # Hz

#     def get_l(self): # luminosity  erg/s
#
#         nu = physics.constants.c / (self.lam * 1E-10)
#
#         return self.Lnu * nu

    def get_Lnu(self, F): # broad band luminosity/erg/s/Hz

        self.Lnu = {f: np.trapz(self.lnu * F[f].T, self.lam) / np.trapz(F[f].T, self.lam) for f in F['filters']}


    def return_Lnu(self, F): # broad band luminosity/erg/s/Hz

        return {f: np.trapz(self.lnu * F[f].T, self.lam) / np.trapz(F[f].T, self.lam) for f in F['filters']}


    def return_Lnu_Window(self, window): # return broad band luminosity/erg/s/Hz in a window defined by window [l_min, l_max]

        s = (self.lam>window[0])&(self.lam<window[1])

        return np.trapz(self.lnu[s], self.lam[s]) / np.trapz(np.ones(len(self.lam[s])), self.lam[s])

    def return_Lnu_lam(self, lam): # return broad band luminosity/erg/s/Hz in a window defined by window [l_min, l_max]

        return np.interp(lam, self.lam, self.lnu)



    def get_fnu(self, cosmo, z, include_IGM = True): # flux nJy, depends on redshift and cosmology

        self.lamz = self.lam * (1. + z)

        self.fnu = 1E23 * 1E9 * self.lnu * (1.+z) / (4 * np.pi * cosmo.luminosity_distance(z).to('cm').value**2) # nJy

        if include_IGM:
            self.fnu *= IGM.madau(self.lamz, z)

    def get_Fnu(self, F): # broad band flux/nJy

        self.Fnu = {f: np.trapz(self.fnu * F[f].T, self.lamz) / np.trapz(F[f].T, self.lamz) for f in F['filters']}

        self.Fnu_array = np.array([self.Fnu[f] for f in F['filters']])

    def return_Fnu(self, F): # broad band flux/nJy

        return {f: np.trapz(self.fnu * F[f].T, self.lamz) / np.trapz(F[f].T, self.lamz) for f in F['filters']}


    # def return_log10Q(self):
    #     """
    #     :return:
    #     """
    #
    #     xi = simps(np.flip(self.lnu[s]/conv), nu(np.flip(lam[s])))/np.interp(nu(1500), nu(lam), fnu)
    #
    #     return np.log10(Q)
    #


    def return_log10Q(self):
        """
        :return:
        """

        llam = self.lnu * c / (self.lam**2*1E-10) # erg s^-1 \AA^-1
        nlam = (llam*self.lam*1E-10)/(h*c) # s^-1 \AA^-1
        s = ((self.lam >= 0) & (self.lam < 912)).nonzero()[0]
        Q = scipy.integrate.simps(nlam[s], self.lam[s])

        return np.log10(Q)

    #
    # def return_log10Q3(self):
    #     """
    #     :return:
    #     """
    #     s = ((self.lam >= 0) & (self.lam < 912)).nonzero()[0]
    #     conv = 1.98644586E-08/(self.lam[s]**2*1E10)
    #     #conv = ((constants.h * constants.c / ((lam[s] * units.AA).to(units.m))).to(units.erg)).value  # alternative using astropy.units and astropy.constants
    #
    #     Q = scipy.integrate.simps(self.lnu[s]/conv, self.lam[s])
    #
    #     return np.log10(Q)




def rebin(l, f, n): # rebin SED [currently destroys original]

    n_len = int(np.floor(len(l)/n))
    l = l[:n_len*n]
    f = f[:n_len*n]
    nl = np.mean(l.reshape(n_len,n), axis=1)
    nf = np.sum(f.reshape(n_len,n), axis=1)/n

    return nl, nf
