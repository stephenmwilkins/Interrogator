

# Create a model SED


import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters



component_plt = False
fesc_plt = False
dust_plt = True
flux_plt = True


# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

# SPS = models.SPS('P2/ModSalpeter_100')
SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')


# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)

sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -2., 'log10M*': 8.})

print('star formation rate: {0}'.format(sfr))





# -------------------------------------------------
# --- plot the stellar, nebular, and total SEDs

if component_plt:

    SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)

    plt.plot(np.log10(SED.stellar.lam), np.log10(SED.stellar.lnu), label = 'stellar')
    plt.plot(np.log10(SED.nebular.lam), np.log10(SED.nebular.lnu), label = 'nebular')
    plt.plot(np.log10(SED.total.lam), np.log10(SED.total.lnu), label = 'total')
    plt.xlim([2.7, 4.])

    mx = np.max(np.log10(SED.total.lnu))
    plt.legend()
    plt.ylim([mx-4., mx+0.3])
    plt.show()






# -------------------------------------------------
# --- plot rest-frame SED for escape fraction of 0.0 and 1.0

if fesc_plt:

    for fesc in [1.0, 0.0]:

        SED = SPS.get_Lnu(sfzh, {'fesc': fesc}, dust = False)

        # plt.plot(np.log10(SED.stellar.lam), np.log10(SED.stellar.lnu))

        plt.plot(np.log10(SED.total.lam), np.log10(SED.total.lnu))

        # -- integrate

        # s = (SED.total.lam>912.)&(SED.total.lam<10000.)
        s = (SED.total.lam<10000.)

        y = SED.total.lnu[s]
        x = SED.total.nu[s]

        Ltot = np.log10(np.trapz(y[::-1], x[::-1]))

        print('fesc: {0} L_tot: {1:.2f}'.format(fesc, Ltot))

        # --- get FUV and NUV luminosities

        filters = ['FAKE.FAKE.1500','FAKE.FAKE.2500']
        F = flare.filters.add_filters(filters, new_lam = SED.total.lam)

        L = SED.total.return_Lnu(F)

        print(L)


    plt.xlim([2.7, 4.])

    mx = np.max(np.log10(SED.total.lnu))
    plt.ylim([mx-4., mx+0.3])
    plt.show()



# -------------------------------------------------
# --- plot the sed with dust

if dust_plt:

    SED = SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': -0.3}, dust = ('simple', {'slope':-1}))

    plt.plot(np.log10(SED.stellar.lam), np.log10(SED.stellar.lnu), label = 'stellar')
    plt.plot(np.log10(SED.total_intrinsic.lam), np.log10(SED.total_intrinsic.lnu), label = 'intrinsic', lw=3, c='0.7')
    plt.plot(np.log10(SED.total.lam), np.log10(SED.total.lnu), label = 'observed')
    plt.xlim([2.7, 4.])

    mx = np.max(np.log10(SED.total_intrinsic.lnu))
    plt.legend()
    plt.ylim([mx-4., mx+0.3])
    plt.show()









# -------------------------------------------------
# --- create observed SED


if flux_plt:

    cosmo = flare.default_cosmo()
    z = 9.0

    SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength

    plt.plot(SED.total.lamz, np.log10(SED.total.fnu), zorder = 1) # plot SED

    filters = flare.filters.HST + flare.filters.Spitzer
    F = flare.filters.add_filters(filters, new_lam = SED.total.lamz) # --- NOTE: need to give it the redshifted


    print(SED.total.lamz)


    SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
    for f in filters:
        plt.scatter(F[f].pivwv(), np.log10(SED.total.Fnu[f]), edgecolor = 'k', zorder = 2, label = f)

    plt.xlim([5000.,50000.])

    mx = np.max(np.log10(SED.total.fnu))
    plt.ylim([mx-4., mx+0.3])
    plt.show()


    # -------------------------------------------------
    # --- save the SED

    # from astropy.io import ascii
    # from astropy.table import Table, Column, MaskedColumn
    #
    # data = Table([SED.total.lamz, SED.total.fnu], names=['observed wavelength/AA','flux/nJy'])
    # ascii.write(data, f'LBG_{z}.dat', overwrite=True)
