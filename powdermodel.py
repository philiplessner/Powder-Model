# coding: utf-8
from __future__ import print_function, division, unicode_literals
from itertools import chain
from IPython.display import HTML
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from utils import format_HTML


class PowderModel(object):

    '''
    Model an anode made of cylinders and pores
    Stores the results of the model as a pandas DataFrame in the instance
    '''

    def __init__(self, anode_props, met, mat_props):
        '''
        Parameters
            a: radius of cylinder after anodization (numpy array)
            met: identity of metal (string)
            Ds: sintered density of anode (g/cc)
            pore_props: {'ratio_ptc': ratio of pores to cylinders,
                         'pore_K: ,
                         'lam': }
        '''
        self.anode_props = anode_props
        self.met = met
        self.mat_props = mat_props

    def _model(self, Vf):
        # Permitivity of free space
        e0 = 8.85e-12
        # Material constants
        X = self.mat_props['Moles metal/Mole oxide']
        K = self.mat_props['dielectric constant']
        rho_m = self.mat_props['density of metal']
        alpha = self.mat_props['microns per volt']
        gamma = self.mat_props['Gamma']
        b0 = self.mat_props['b0']
        # Anode properties
        a = self.anode_props['a']
        Ds = self.anode_props['Ds']
        ratio_ptc = self.anode_props['ratio_ptc']
        pore_K = self.anode_props['pore_K']
        lam = self.anode_props['lam']
        # Cylinder Model
        t_oxide = b0 + alpha * Vf
        b = a + t_oxide
        Rsq = X * gamma * (b ** 2 - a ** 2) + a ** 2
        R = np.sqrt(Rsq)
        Dnm = 2. * R * 1000.
        CV_gcyl = 2. * e0 * K * Vf / (Rsq * rho_m * np.log(b / a))
        rho_corr_factor = (0.8529 - Ds / rho_m * 0.712 * 1.e-12)
        corrCV_gcyl = CV_gcyl * rho_corr_factor
        # Cylinder + pore model
        Rp = R * ratio_ptc
        Pnm = 2. * Rp * 1000.
        Rpsq = Rp ** 2
        ap = (-2. * (1. - 2 * gamma) * t_oxide + np.sqrt(4. * ((1. - 2. * gamma) ** 2) * (t_oxide ** 2)
                                                         - 4. * ((1 - 2. * gamma) * (t_oxide ** 2) - Rpsq))) / 2.
        bp = ap + t_oxide
        CV_gpore_cyl = (corrCV_gcyl * Rsq) / (Rsq + (pore_K - 1.) * Rpsq * lam)
        # CV/cc (if lam == 0, defaults to cylinder model)
        CV_cc = Ds * CV_gpore_cyl
        return {'a': a,
                'Vf': Vf,
                'Ds': Ds,
                'b': b,
                'Rsq': Rsq,
                'R': R,
                'Dnm': Dnm,
                'CV per g cyl': CV_gcyl,
                'density corr factor': rho_corr_factor,
                'corrCV per g cyl': corrCV_gcyl,
                'Rp': Rp,
                'Pnm': Pnm,
                'Rpsq': Rpsq,
                'ap': ap,
                'bp': bp,
                'CV per g por_cyl': CV_gpore_cyl,
                'CV per cc': CV_cc}

    def run_model(self, V):
        self.df_overall = reduce(pd.DataFrame.append,
                                 [pd.DataFrame(self._model(Vf)) for Vf in V],
                                 pd.DataFrame())

    def get_modeldata(self):
        return self.df_overall

    def get_HTMLtable(self):
        return HTML(self.df_overall.to_html(index=False,
                                            **format_HTML([['a', '{:.6f}'],
                                                           ['Vf', '{:.0f}'],
                                                           ['Ds', '{:.2f}'],
                                                           ['Dnm', '{:.1f}'],
                                                           ['CV per g cyl',
                                                               '{:,.0f}'],
                                                           ['density corr factor',
                                                               '{:.4f}'],
                                                           ['corrCV per g cyl',
                                                               '{:,.0f}'],
                                                           ['CV per cc',
                                                            '{:,.0f}']])))

    def make_plots(self):
        V = np.unique(self.df_overall['Vf'])
        gb_vf = self.df_overall.groupby('Vf')
        majorFormatter = FuncFormatter(lambda x, pos: '{:,.0f}'.format(x))
        fontsize = 18
        fig, ax = plt.subplots(2, 2, figsize=(20, 9))
        for name, group in gb_vf:
            ax[0, 0].plot(group['Dnm'],
                          group['corrCV per g cyl'],
                          label='$V_f=$' + unicode(name))
            ax[1, 0].plot(group['Dnm'],
                          group['CV per cc'],
                          label='$V_f=$' + unicode(name))
        ax[0, 0].set_ylabel('CV/g', fontsize=fontsize)
        ax[0, 0].set_xlabel('Cylinder Diameter (nm)', fontsize=fontsize)
        ax[0, 0].legend()
        fig.suptitle(' '.join([self.met,
                               'Anode Density',
                               unicode(self.anode_props['Ds'])]),
                     fontsize=fontsize)
        ax[1, 0].set_xlabel('Cylinder Diameter (nm)', fontsize=fontsize)
        ax[1, 0].set_ylabel('CV/cc', fontsize=fontsize)
        ax[1, 0].legend()
        ax[0, 1].scatter(1 / V, gb_vf['corrCV per g cyl'].max())
        CVg_fit = np.ma.polyfit(1. / V, gb_vf['corrCV per g cyl'].max(), 1)
        ax[0, 1].plot(1. / V, np.polyval(CVg_fit, 1. / V))
        ax[0, 1].set_xlabel('1/$V_f$ ($V^{-1}$)', fontsize=fontsize)
        ax[0, 1].set_ylabel('$CV/g_{max}$', fontsize=fontsize)
        ax[0, 1].text(0.3 * ax[0, 1].get_xlim()[1],
                      0.9 * ax[0, 1].get_ylim()[1],
                      '(CV/g)max = {:,.0f}/Vf + {:,.0f}'.format(CVg_fit[0],
                                                                CVg_fit[1]))
        ax[1, 1].scatter(1 / V, gb_vf['CV per cc'].max())
        CVcc_fit = np.ma.polyfit(1. / V, gb_vf['CV per cc'].max(), 1)
        ax[1, 1].plot(1. / V, np.polyval(CVcc_fit, 1. / V))
        ax[1, 1].text(0.3 * ax[1, 1].get_xlim()[1],
                      0.9 * ax[1, 1].get_ylim()[1],
                      '(CV/cc)max = {:,.0f}/Vf + {:,.0f}'.format(CVcc_fit[0],
                                                                 CVcc_fit[1]))
        ax[1, 1].set_xlabel('1/$V_f$ ($V^{-1}$)', fontsize=fontsize)
        ax[1, 1].set_ylabel('$CV/cc_{max}$', fontsize=fontsize)
        for e in chain.from_iterable(ax):
            e.yaxis.set_major_formatter(majorFormatter)
            plt.setp(e.get_xticklabels(), fontsize=fontsize - 4)
            plt.setp(e.get_yticklabels(), fontsize=fontsize - 4)
        plt.show()
        return fig, ax


def compare_metals(metals_data, labels, Vf):
    '''
    Compare CV/g and CV/cc for different metals at a given formation voltage
    Parameters
        metals_data: pandas groups at voltage
        labels: labels for legend
        Vf: voltage corresponding to the groups voltage
    '''
    majorFormatter = FuncFormatter(lambda x, pos: '{:,.0f}'.format(x))
    fontsize = 18
    fig, ax = plt.subplots(1, 2, figsize=(20, 5))
    fig.suptitle(''.join(['$V_f=$',  unicode(Vf), '$V$']), fontsize=fontsize)
    for metal_data, label in zip(metals_data, labels):
        ax[0].plot(
            metal_data['Dnm'], metal_data['corrCV per g cyl'], label=label)
        ax[1].plot(metal_data['Dnm'], metal_data['CV per cc'], label=label)
    for e in ax:
        e.set_xlabel('Cylinder Diameter (nm)', fontsize=fontsize)
        e.set_xlim(left=0.)
        e.yaxis.set_major_formatter(majorFormatter)
        e.legend(loc='best')
        plt.setp(e.get_xticklabels(), fontsize=fontsize - 4)
        plt.setp(e.get_yticklabels(), fontsize=fontsize - 4)
    ax[0].set_ylabel('CV/g', fontsize=fontsize)
    ax[1].set_ylabel('CV/cc', fontsize=fontsize)
    plt.show()
    return fig, ax
