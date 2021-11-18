
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import flare
import flare.filters
import flare.observatories as obs

from interrogator.sed import models
import interrogator.sed.sfzh


z = 8

cosmo = flare.default_cosmo()

SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')
filters = obs.Hubble.WFC3.f() + obs.Spitzer.IRAC.f()
print(filters)
F = flare.filters.add_filters(filters, new_lam = SPS.lam * (1. + z)) # --- NOTE: need to give it the redshifted

cxf = ['Hubble.WFC3.f125w', 'Spitzer.IRAC.ch1']
cyf = ['Spitzer.IRAC.ch1', 'Spitzer.IRAC.ch2']



# ---- constant

for log10tau_V in [-2, -1, -0.5, 0.0]:

    print(f'------- log10tau_V: {log10tau_V}')

    for a0 in [7.,8.,9.]:

        sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'], {'log10_duration': a0, 'log10Z': -2., 'log10M*': 8.})

        SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))

        # --- create observed SED

        SED.total.get_fnu(cosmo, z=8.0) # calculate observer frame wavelength

        SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)

        fluxes = SED.total.Fnu

        cx = -2.5*np.log10(fluxes[cxf[0]]/fluxes[cxf[1]])
        cy = -2.5*np.log10(fluxes[cyf[0]]/fluxes[cyf[1]])

        print(a0,cx,cy, SED.A1500())
