# coding: utf-8
from __future__ import print_function, division, unicode_literals
import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from toolz.dicttoolz import merge_with


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
                'Vf': np.repeat(Vf, len(a)),
                'Ds': np.repeat(Ds, len(a)),
                'b': b,
                'Rsq': Rsq,
                'R': R,
                'Dnm': Dnm,
                'CV per g cyl': CV_gcyl,
                'density corr factor': np.repeat(rho_corr_factor, len(a)),
                'corrCV per g cyl': corrCV_gcyl,
                'Rp': Rp,
                'Pnm': Pnm,
                'Rpsq': Rpsq,
                'ap': ap,
                'bp': bp,
                'CV per g por_cyl': CV_gpore_cyl,
                'CV per cc': CV_cc}

    def run_model(self, V):
        for i, Vf in enumerate(reversed(V)):
            cap = self._model(Vf) if i == 0 else merge_with(np.concatenate,
                                                            self._model(Vf),
                                                            cap)
        return cap


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
