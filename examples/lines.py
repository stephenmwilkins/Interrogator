




import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh




# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

# SPS = models.SPS('P2/ModSalpeter_100')
line_modeller = models.lines('BPASSv2.2.1.binary/ModSalpeter_300')


print(list(line_modeller.grid.keys()))
print(list(line_modeller.grid['HI4861'].keys()))


# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)

sfzh, sfr = interrogator.sed.sfzh.constant(line_modeller.grid['log10age'], line_modeller.grid['log10Z'] , {'log10_duration': 7., 'log10Z': -2., 'log10M*': 8.})

print('star formation rate: {0}'.format(sfr))

fesc = 0.5
Hbeta = line_modeller.get_info(sfzh, 'HI4861', {'fesc': fesc}, dust_model = ('just_gas'))

print(Hbeta)


# -------------------------------------------------
# --- all lines




all_line_info = line_modeller.get_all_info_cols(sfzh, {'fesc': fesc}, dust_model = ('just_gas'))

print(all_line_info)
