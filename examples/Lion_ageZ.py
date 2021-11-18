
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from interrogator.sed import models
import interrogator.sed.sfzh



make_nice_plot = True



SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')

log10_durations =  SPS.grid['log10age']




def get_log10xi_ion(log10_durations, log10Z, verbose = False):

    log10xi_ion = np.zeros(len(log10_durations))
    print(log10Z)

    for i, log10_duration in enumerate(log10_durations):

        sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'], {'log10_duration': log10_duration, 'log10Z': log10Z, 'log10M*': 8.})


        # --- get L1500 luminosity
        SED = SPS.get_Lnu(sfzh, {'fesc': 1.0})
        # log10L1500 = np.log10(SED.total.return_Lnu_lam(1500.)) # --- using a single wavelength
        log10L1500 = np.log10(SED.total.return_Lnu_Window([1400.,1600.])) # --- using a spectral window

        # --- get the ionising photon "luminosity"
        log10Q = SPS.get_log10Q(sfzh)

        # --- calculate the ionising photon production efficency
        log10xi_ion[i] = log10Q - log10L1500

        if verbose: print(log10_duration, log10xi_ion)

    return log10xi_ion









if make_nice_plot:

    import flare.plt as fplt
    import cmasher as cmr
    import matplotlib as mpl
    import matplotlib.cm as cm

    cmap = cmr.sapphire

    # log10Z_min = SPS.grid['log10Z'][0]
    # log10Z_max = SPS.grid['log10Z'][-1]
    # log10Z_ref = -2.
    #
    # log10xi_ion_min = get_log10xi_ion(log10_durations, log10Z_min)
    # log10xi_ion_max = get_log10xi_ion(log10_durations, log10Z_max)
    # log10xi_ion_ref = get_log10xi_ion(log10_durations, log10Z_ref)
    #
    # fig, ax = fplt.simple()
    #
    # ax.fill_between(log10_durations, log10xi_ion_min, log10xi_ion_max, color='k', alpha=0.1, label = rf'$\rm log_{{10}}(Z)=[{log10Z_min:.2f},{log10Z_max:.2f}] $')
    #
    # ax.plot(log10_durations, log10xi_ion_ref, c='k', label = rf'$\rm log_{{10}}(Z)={log10Z_ref:.2f} $')
    #
    # ax.legend(fontsize=8)
    #
    # ax.set_ylabel(fplt.fancy(r'log{10}(\xi_{ion}/s^{-1}/erg s^{-1} Hz^{-1})'))
    # ax.set_xlabel(fplt.fancy(r'log_{10}(SF duration/Myr)'))
    #
    # ax.set_xlim([log10_durations[0],log10_durations[-1]])
    # # ax.set_ylim([-0.499, 3.99])
    #
    # ax.legend(fontsize = 8)
    #
    # fig.savefig(f'Lion_ageZ.pdf')






    fig, ax, cax = fplt.simple_wcbar()


    norm = mpl.colors.Normalize(vmin=SPS.grid['log10Z'][0], vmax=SPS.grid['log10Z'][-1])

    for log10Z in SPS.grid['log10Z']:

        c = cmap(norm(log10Z))

        log10xi_ion = get_log10xi_ion(log10_durations, log10Z)

        # ax.plot(log10_durations, log10xi_ion, c=c, alpha = 0.8, label = rf'$\rm {log10Z:.2f} $')
        ax.plot(log10_durations, log10xi_ion, c=c, alpha = 0.8)




    ax.legend(fontsize=8)

    ax.set_ylabel(fplt.fancy(r'log{10}(\xi_{ion}/s^{-1}/erg s^{-1} Hz^{-1})'))
    ax.set_xlabel(fplt.fancy(r'log_{10}(SF duration/Myr)'))

    ax.set_xlim([log10_durations[0],log10_durations[-1]])
    # ax.set_ylim([-0.499, 3.99])

    # ax.legend(fontsize = 8)


    cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    cmapper.set_array([])
    cbar = fig.colorbar(cmapper, cax=cax, orientation='vertical')
    cbar.set_label(r'$\rm log_{10}Z $')

    fig.savefig(f'Lion_ageZ.pdf')



else:

    log10Z_min = SPS.grid['log10Z'][0]
    log10Z_max = SPS.grid['log10Z'][-1]
    log10Z_ref = -2.

    log10xi_ion_min = get_log10xi_ion(log10_durations, log10Z_min)
    log10xi_ion_max = get_log10xi_ion(log10_durations, log10Z_max)
    log10xi_ion_ref = get_log10xi_ion(log10_durations, log10Z_ref)

    plt.plot(log10_durations, log10xi_ion_min, label = rf"$\rm log10(Z)={log10Z_min:.2f} $")
    plt.plot(log10_durations, log10xi_ion_ref, label = rf"$\rm log10(Z)={log10Z_ref:.2f} $")
    plt.plot(log10_durations, log10xi_ion_max, label = rf"$\rm log10(Z)={log10Z_max:.2f} $")
    plt.legend()
    plt.show()
