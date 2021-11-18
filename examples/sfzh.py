
# Produce different star formation and metal enrichment histories


import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters


# SPS = models.SPS('P2/ModSalpeter_100')
SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')


BW = interrogator.sed.sfzh.bin_width(SPS.grid['log10age'], SPS.grid['log10Z'])

sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -2., 'log10M*': 8.})


plt.imshow(sfzh)
plt.show()

plt.imshow(sfzh/BW)
plt.show()


sfzh_q, sfr = interrogator.sed.sfzh.quenched(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_quenched_duration': 7.,'log10_total_duration': 8., 'log10Z': -2., 'log10M*': 8.})

plt.imshow(sfzh_q)
plt.show()



sfzh_exp, sfr = interrogator.sed.sfzh.exponential(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 9.,'tau_Myr': 100., 'log10Z': -2., 'log10M*': 8.})

plt.imshow(sfzh_exp)
plt.show()
plt.imshow(sfzh_exp/BW)
plt.show()
