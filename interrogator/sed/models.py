

import numpy as np
import pickle

from . import core
# from ..core import *

import flare

import copy

from . import dust_curves

class empty: pass



def beta(lam, slope, normalisation, normalisation_wavelength = 1500., include_ISM = True):

    model = core.sed(lam)

    model.lnu = normalisation*(lam/normalisation_wavelength)**(slope + 2.0)

    if include_ISM: model.lnu[lam<912.] = 0.0 # add ISM absorption

    return model



class SED():

    def A(self,l):

        return -2.5*np.log10(np.interp(l, self.total.lam, self.total.lnu)/np.interp(l, self.total_intrinsic.lam, self.total_intrinsic.lnu))

    def A1500(self):

        return self.A(1500.)




class SPS():

    def __init__(self, grid, path_to_SPS_grid = '/data/SPS/nebular/3.0/'):

        self.grid = pickle.load(open(flare.FLARE_dir + path_to_SPS_grid + grid + '/nebular.p','rb'), encoding='latin1')

        self.lam = self.grid['lam']


    def get_Lnu(self, sfzh, SED_p, dust = False, fast = False):

        sed = SED()


        SFZH = np.expand_dims(sfzh, axis=2)

        # --- calculate intrinsic SED

        sed.stellar = core.sed(self.lam)
        sed.stellar.lnu = np.sum(self.grid['stellar'] * SFZH, axis=(0,1))
        sed.stellar.lnu[self.lam<912.] *= SED_p['fesc']

        sed.nebular = core.sed(self.lam)
        sed.nebular.lnu = np.sum(self.grid['nebular'] * SFZH, axis=(0,1))

        sed.total = core.sed(self.lam)
        sed.total.lnu = sed.stellar.lnu + (1.-SED_p['fesc'])*sed.nebular.lnu


        if dust:

            # --- intrinsic SEDs
            sed.stellar_intrinsic = copy.deepcopy(sed.stellar)
            sed.nebular_intrinsic = copy.deepcopy(sed.nebular)
            sed.total_intrinsic = copy.deepcopy(sed.total)


            dust_model, dust_model_params = dust


            # --- create SED of young and old components

            if dust_model == 'pacman':

                dmp = dust_model_params

                # --- see Ciaran's thesis example for how to use this

                i = np.where(self.grid['log10age']>dmp['log10age_BC'])[0][0]

                sfzh_young = copy.copy(sfzh)
                sfzh_young[i:,:] = 0.0

                sfzh_old = copy.copy(sfzh)
                sfzh_old[:i,:] = 0.0

                SFZH_young = np.expand_dims(sfzh_young, axis=2)
                SFZH_old = np.expand_dims(sfzh_old, axis=2)


                tau_BC = 10**(dmp['tau_ISM_to_BC']*SED_p['log10tau_V']) * dust_curves.simple(params = {'slope': dmp['alpha_BC']}).tau(self.lam)

                tau_ISM = 10**(SED_p['log10tau_V']) * dust_curves.simple(params = {'slope': dmp['alpha_ISM']}).tau(self.lam)

                T_young = np.exp(-tau_BC)*np.exp(-tau_ISM)
                T_old = np.exp(-tau_ISM)


                # --- young stellar light which escapes with no gas/dust reprocessing
                sed.stellar_young_esc = core.sed(self.lam)
                sed.stellar_young_esc.lnu = SED_p['fesc'] * np.sum(self.grid['stellar'] * SFZH_young, axis=(0,1))

                sed.stellar_young_nesc = core.sed(self.lam)
                sed.stellar_young_nesc.lnu = T_young*(1-SED_p['fesc']) * np.sum(self.grid['stellar'] * SFZH_young, axis=(0,1))
                sed.stellar_young_nesc.lnu[self.lam<912.] *= SED_p['fesc']


                sed.nebular_young_nesc = core.sed(self.lam)
                sed.nebular_young_nesc.lnu = T_young*(1-SED_p['fesc']) * np.sum(self.grid['nebular'] * SFZH_young, axis=(0,1))

                sed.total_young_esc = core.sed(self.lam)
                sed.total_young_esc.lnu = sed.stellar_young_esc.lnu # should be same

                sed.total_young = core.sed(self.lam)
                sed.total_young.lnu = sed.total_young_esc.lnu + sed.stellar_young_nesc.lnu + sed.nebular_young_nesc.lnu

                # --- old stellar light which escapes with no gas/dust reprocessing
                sed.stellar_old_esc = core.sed(self.lam)
                sed.stellar_old_esc.lnu = SED_p['fesc'] * np.sum(self.grid['stellar'] * SFZH_old, axis=(0,1))

                sed.stellar_old_nesc = core.sed(self.lam)
                sed.stellar_old_nesc.lnu = T_old*(1-SED_p['fesc']) * np.sum(self.grid['stellar'] * SFZH_old, axis=(0,1))
                sed.stellar_old_nesc.lnu[self.lam<912.] *= SED_p['fesc']

                sed.nebular_old_nesc = core.sed(self.lam)
                sed.nebular_old_nesc.lnu = T_old*(1-SED_p['fesc']) * np.sum(self.grid['nebular'] * SFZH_old, axis=(0,1))

                sed.total_old_esc = core.sed(self.lam)
                sed.total_old_esc.lnu = sed.stellar_young_esc.lnu # should be same

                sed.total_old = core.sed(self.lam)
                sed.total_old.lnu = sed.total_old_esc.lnu + sed.stellar_old_nesc.lnu + sed.nebular_old_nesc.lnu

                # --- total
                sed.total = core.sed(self.lam)
                sed.total.lnu = sed.total_young.lnu + sed.total_old.lnu


            elif dust_model == 'CF00':

                print('not yet implemented')

            elif dust_model == 'simple':

                tau = 10**(SED_p['log10tau_V']) * getattr(dust_curves, dust_model)(params = dust_model_params).tau(self.lam)

                T = np.exp(-tau)

                sed.stellar.lnu *= T
                sed.nebular.lnu *= T
                sed.total.lnu *= T


            else:

                print('dust model not implemented')

        return sed



    def get_Q(self, SFZH):
        return np.log10(np.sum(10**self.grid['log10Q'] * SFZH, axis=(0,1)))

    def get_log10Q(self, SFZH):
        return self.get_Q(SFZH)





# --- this needs updating to do lines, will also need dust

class lines():

    def __init__(self, grid, path_to_SPS_grid = flare.FLARE_dir + '/data/SPS/nebular/3.0/'):

        self.grid = pickle.load(open(path_to_SPS_grid + grid + '/lines.p','rb'), encoding='latin1')



    def get_all_info(self, SFZH, SED_p, dust = False):

        return {line: self.get_info(SFZH, line, SED_p, dust = dust) for line in self.grid['lines']}


    # --- return above as columns
    def get_all_info_cols(self, SFZH, SED_p, dust = False):

        line_info = self.get_all_info(SFZH, SED_p, dust = dust)
        print(line_info.keys())

        line_info_col = {k:[] for k in ['line', 'lam', 'EW', 'luminosity']}

        for line in self.grid['lines']:
            line_info_col['line'].append(line)
            for col in ['lam', 'EW', 'luminosity']:
                line_info_col[col].append(line_info[line][col])

        for col in ['line','lam', 'EW', 'luminosity']:
            line_info_col[col] = np.array(line_info_col[col])

        return line_info_col


    def get_info(self, SFZH, line, SED_p, dust = False):  # not this was changed to  match the usage in SED modeller, this may have broken something.


        if dust:
            if len(dust) == 2:
                dust_model, dust_model_params = dust
            else:
                dust_model = dust
        else:
            dust_model = False


        line_info = {}

        # --- wavelength

        line_info['lam'] = lam = self.grid[line]['lam'] #AA

        # --- these assume fesc = 0.0
        line_info['max_luminosity'] = np.sum(10**self.grid[line]['luminosity'] * SFZH, axis=(0,1))
        line_info['max_total_continuum'] = np.sum(self.grid[line]['total_continuum'] * SFZH, axis=(0,1))
        line_info['max_EW'] = line_info['max_luminosity']/(3E8*line_info['max_total_continuum']/((line_info['lam'])**2*1E-10))

        if dust_model:

            if dust_model in ['just_gas', 'simple', 'pacman']:

                line_info['intrinsic_luminosity'] = (1-SED_p['fesc'])*np.sum(10**self.grid[line]['luminosity'] * SFZH, axis=(0,1))
                line_info['intrinsic_nebular_continuum'] = (1-SED_p['fesc'])*np.sum(self.grid[line]['nebular_continuum'] * SFZH, axis=(0,1))
                line_info['intrinsic_stellar_continuum'] = np.sum(self.grid[line]['stellar_continuum'] * SFZH, axis=(0,1))
                line_info['intrinsic_total_continuum'] =  line_info['intrinsic_nebular_continuum'] + line_info['intrinsic_stellar_continuum']
                line_info['intrinsic_EW'] = line_info['intrinsic_luminosity']/(3E8*line_info['intrinsic_total_continuum']/((line_info['lam'])**2*1E-10))

                if dust_model == 'just_gas':

                    line_info['luminosity'] = copy.copy(line_info['intrinsic_luminosity'])
                    line_info['total_continuum'] = copy.copy(line_info['intrinsic_total_continuum'])
                    line_info['EW'] = copy.copy(line_info['intrinsic_EW'])


            if dust_model == 'simple':

                tau = 10**(SED_p['log10tau_V']) * getattr(dust_curves, dust_model)(params = dust_model_params).tau(lam)
                T = np.exp(-tau)

                line_info['attenuated_luminosity'] = T * line_info['intrinsic_luminosity']
                line_info['attenuated_nebular_continuum'] =  T * line_info['intrinsic_nebular_continuum']
                line_info['attenuated_stellar_continuum'] = T * line_info['intrinsic_stellar_continuum']
                line_info['attenuated_total_continuum'] =  line_info['attenuated_nebular_continuum'] + line_info['attenuated_stellar_continuum']
                line_info['attenuated_EW'] = line_info['attenuated_luminosity']/(3E8*line_info['attenuated_total_continuum']/((line_info['lam'])**2*1E-10))

                line_info['luminosity'] = line_info['attenuated_luminosity']
                line_info['nebular_continuum'] = line_info['attenuated_nebular_continuum']
                line_info['stellar_continuum'] = line_info['attenuated_stellar_continuum']
                line_info['total_continuum'] =  line_info['attenuated_total_continuum']
                line_info['EW'] = line_info['attenuated_EW']



            if dust_model == 'pacman':

                dmp = dust_model_params

                i = np.where(self.grid['log10age']>dmp['log10age_BC'])[0][0]

                sfzh_young = copy.copy(sfzh)
                sfzh_young[i:,:] = 0.0

                sfzh_old = copy.copy(sfzh)
                sfzh_old[:i,:] = 0.0

                SFZH_young = np.expand_dims(sfzh_young, axis=2)
                SFZH_old = np.expand_dims(sfzh_old, axis=2)

                line_info['intrinsic_young_luminosity'] = (1-SED_p['fesc'])*np.sum(10**self.grid[line]['luminosity'] * SFZH_young, axis=(0,1))
                line_info['intrinsic_young_nebular_continuum'] = (1-SED_p['fesc'])*np.sum(self.grid[line]['nebular_continuum'] * SFZH_young, axis=(0,1))
                line_info['intrinsic_young_stellar_continuum'] = np.sum(self.grid[line]['stellar_continuum'] * SFZH_young, axis=(0,1))
                line_info['intrinsic_young_total_continuum'] =  line_info['intrinsic_young_nebular_continuum'] + line_info['intrinsic_young_stellar_continuum']

                line_info['intrinsic_old_luminosity'] = (1-SED_p['fesc'])*np.sum(10**self.grid[line]['luminosity'] * SFZH_old, axis=(0,1))
                line_info['intrinsic_old_nebular_continuum'] = (1-SED_p['fesc'])*np.sum(self.grid[line]['nebular_continuum'] * SFZH_old, axis=(0,1))
                line_info['intrinsic_old_stellar_continuum'] = np.sum(self.grid[line]['stellar_continuum'] * SFZH_old, axis=(0,1))
                line_info['intrinsic_old_total_continuum'] =  line_info['intrinsic_old_nebular_continuum'] + line_info['intrinsic_old_stellar_continuum']

                tau_BC = 10**(dmp['tau_ISM_to_BC']*SED_p['log10tau_V']) * dust_curves.simple(params = {'slope': dmp['alpha_BC']}).tau(self.lam)
                tau_ISM = 10**(SED_p['log10tau_V']) * dust_curves.simple(params = {'slope': dmp['alpha_ISM']}).tau(self.lam)

                T_young = np.exp(-tau_BC)*np.exp(-tau_ISM)
                T_old = np.exp(-tau_ISM)

                line_info['attenuated_young_luminosity'] = T_young * line_info['intrinsic_young_luminosity']
                line_info['attenuated_young_nebular_continuum'] = T_young * line_info['intrinsic_young_nebular_continuum']
                line_info['attenuated_young_stellar_continuum'] = T_young * line_info['intrinsic_young_stellar_continuum']
                line_info['attenuated_young_total_continuum'] = T_young * line_info['intrinsic_young_total_continuum']

                line_info['attenuated_old_luminosity'] = T_young * line_info['intrinsic_old_luminosity']
                line_info['attenuated_old_nebular_continuum'] = T_young * line_info['intrinsic_old_nebular_continuum']
                line_info['attenuated_old_stellar_continuum'] = T_young * line_info['intrinsic_old_stellar_continuum']
                line_info['attenuated_old_total_continuum'] = T_young * line_info['intrinsic_old_total_continuum']

                line_info['attenuated_luminosity'] = line_info['attenuated_young_luminosity'] + line_info['attenuated_old_luminosity']
                line_info['attenuated_nebular_continuum']= line_info['attenuated_young_nebular_continuum'] + line_info['attenuated_old_nebular_continuum']
                line_info['attenuated_stellar_continuum'] = line_info['attenuated_young_stellar_continuum']+ line_info['attenuated_old_stellar_continuum']
                line_info['attenuated_total_continuum'] = line_info['attenuated_young_total_continuum'] + line_info['attenuated_old_total_continuum']
                line_info['attenuated_EW'] = line_info['attenuated_luminosity']/(3E8*line_info['attenuated_total_continuum']/((line_info['lam'])**2*1E-10))

                line_info['luminosity'] = copy.copy(line_info['luminosity'])
                line_info['nebular_continuum'] = copy.copy(line_info['nebular_continuum'])
                line_info['stellar_continuum'] = copy.copy(line_info['stellar_continuum'])
                line_info['total_continuum'] = copy.copy(line_info['total_continuum'])
                line_info['EW'] = copy.copy(line_info['attenuated_EW'])

        return line_info
