# coding: utf-8
from __future__ import print_function, division, unicode_literals
import csv
import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from operator import add, mul
from toolz import thread_first
from toolz.dicttoolz import merge_with


class PowderModel(object):

    '''
    Model an anode made of cylinders and pores
    '''

    def __init__(self, anode_props, met, mat_props):
        '''
        Parameters
            anode_props: dictionary
                a: radius of cylinders after anodization (numpy array)
                Ds: density of sintered anode (g/cc)
                ratio_ptc: ratio of pores to cylinders
                pore_K:
                lam:
            met: identity of metal (string)
            mat_props: dictionary
                dielectric constant: relative dielectric constant
                density of oxide:
                density of metal:
                microns per volt: formation constant
                b0: thickness of native oxide film (microns)
                Molecular Weight of Metal:
                Molecular Weight oxide:
                Gamma:
                X: Moles of metal per mole of metal oxide
                Pilling:
        '''
        self.anode_props = anode_props
        self.met = met
        self.mat_props = mat_props
        # cap holds the results of the model after run_model
        self.cap = None

    def _model(self, Vf):
        '''
        Calculate capacitor properties at a given formation voltage
        Parameter
            Vf: formation voltage (volts)
        Returns
            dictionary containing input and calculated properties
                a: (input property) inner radius of cylinders after formation
                Vf: (input property) formation voltage (volts)
                Ds: (input property) sintered density of anode (g/cc)
                b:  outer radius of cylinders after formation (microns)
                Rsq: square of radius of cylinder before formation
                R: radius of cylinder before formation (microns)
                Dnm: diameter of cylinder before formation (nm)
                CV per g cyl: uC/g of cylinder not corr for sintered density
                density corr factor: corr factor loss of surface area
                corrCV per g cyl: uC/g after correction for density
                Rp: radius of pore (microns)
                Pnm: diameter of pore (nm)
                Rpsq: radius of pore squared
                ap:
                bp:
                CV per g pore_cyl: uC/g for cylinder + pore model
                CV per cc: uC/cc
        '''
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
        '''
        Run the model over a range of voltages.
        Results are stored in cap
        Parameter
            V: ndarray of voltages
        '''
        for i, Vf in enumerate(reversed(V)):
            self.cap = self._model(Vf) if i == 0 else merge_with(
                np.concatenate, self._model(Vf), self.cap)

    def make_plots(self):
        if not self.cap:
            print('\nRun model before plotting\n')
        # fontsize=18
        K = self.mat_props['dielectric constant']
        alpha = self.mat_props['microns per volt']
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
        majorFormatter = FuncFormatter(lambda x, pos: '{:,.0f}'.format(x))
        for z in np.unique(self.cap['Vf']):
            x = self.cap['Dnm'][self.cap['Vf'] == z]
            ax[0].plot(x,
                       self.cap['corrCV per g cyl'][self.cap['Vf'] == z],
                       label=' '.join(['$V_f =$', unicode(z), 'V']))
            ax[1].plot(x,
                       self.cap['CV per cc'][self.cap['Vf'] == z],
                       label=' '.join(['$V_f =$', unicode(z), 'V']))
        for e in ax:
            e.yaxis.set_major_formatter(majorFormatter)
            e.set_xlabel('Cylinder Diameter (nm)')
            e.legend()
        ek = '{:.1f}'.format(K / (alpha * 1000.))
        Pilling = '{:.3f}'.format(self.mat_props['Pilling'])
        ax[0].set_title(' '.join([self.met, '(',
                                  '$\epsilon /k=$', ek,
                                  ',', 'Pilling=',
                                  Pilling, ')']))
        ad = '{:.2f}'.format(self.anode_props['Ds'])
        ax[1].set_title(' '.join([self.met, '(', '$D_s=$', ad, ')']))
        ax[0].set_ylabel('CV/g ($\mu C/g$)')
        ax[1].set_ylabel('CV/cc ($\mu C/cc$)')
        plt.tight_layout()
        plt.show()
        plt.clf()

    def get_dataV(self, V):
        '''
        After running model, get results at a given voltage.
        Parameter
            V: voltage (volts)
        Returns
            dictionary whose values for each key are ndarrays
            of model results
        '''
        return {key: self.cap[key][self.cap['Vf'] == V]
                for key in self.cap}

    def v_rolloff(self, D):
        '''
        Compute CV/g rolloff with Voltage for a given
        diameter cylinder
        Parameter
            D: cylinder diameter (nm)
        '''
        Vs = np.unique(self.cap['Vf'])
        y = np.zeros(len(Vs))
        for i, V in enumerate(Vs):
            d = self.get_dataV(V)
            Ds_gt = d['Dnm'][d['Dnm'] > D]
            Ds_lt = d['Dnm'][d['Dnm'] < D]
            # Check if that diameter is available for voltage
            if Ds_gt.any() and Ds_lt.any():
                # Use linear interpolation
                D_gt = Ds_gt[0]
                CVg_gt = d['corrCV per g cyl'][d['Dnm'] > D][0]
                D_lt = Ds_lt[-1]
                CVg_lt = d['corrCV per g cyl'][d['Dnm'] < D][-1]
                m = (CVg_gt - CVg_lt) / (D_gt - D_lt)
                CVg = m * (D - D_lt) + CVg_lt
                # print(V, D, CVg)
                y[i] = CVg
            else:
                y[i] = np.nan
        return {'D': D,
                'Vs': Vs,
                'CVg': y}

    def write_csv(self, fpath):
        '''
        Write model results to a csv file.
        Parameters
            fpath: path to csv file
        '''
        col_order = ['a', 'Vf', 'Ds', 'b', 'Rsq', 'R', 'Dnm', 'CV per g cyl',
                     'density corr factor', 'corrCV per g cyl', 'Rp', 'Pnm',
                     'Rpsq', 'ap', 'bp', 'CV per g por_cyl', 'CV per cc']
        with open(fpath, 'w') as f:
            f.write(','.join([key for key in col_order]))
            f.write('\n')
            for i in xrange(0, self.cap['Ds'].shape[0]):
                s = ','.join([unicode(self.cap[key][i]) for key in col_order])
                f.write(''.join([s, '\n']))


class FlakeModel(object):

    '''
    Model an anode made of flakes
    '''

    def __init__(self, anode_props, met, mat_props):
        '''
        Parameters
            anode_props: dictionary
                r0: radius of flakes before anodization
                t0: thickness of flakes
                L: array of neck thickness
                Ds: density of sintered anode (g/cc)
            met: identity of metal (string)
            mat_props: dictionary
                dielectric constant: relative dielectric constant
                density of oxide:
                density of metal:
                microns per volt: formation constant
                b0: thickness of native oxide film (microns)
                Molecular Weight of Metal:
                Molecular Weight oxide:
                Gamma:
                X: Moles of metal per mole of metal oxide
                Pilling:
        '''
        self.anode_props = anode_props
        self.met = met
        self.mat_props = mat_props
        # cap holds the results of the model after run_model
        self.cap = None

    def _flakemodel(self, Vf):
        e0 = 8.85e-12
        r0 = self.anode_props['r0']
        t0 = self.anode_props['t0']
        Ds = self.anode_props['Ds']
        L = self.anode_props['L']
        K = self.mat_props['dielectric constant']
        rho_m = self.mat_props['density of metal']
        alpha = self.mat_props['microns per volt']
        P = self.mat_props['Pilling']
        # Propeties of sintered anode
        r = thread_first(L ** 2 * t0 * (L - t0) / 2.,
                         (add, np.sqrt(t0 ** 3 * (4 * L ** 3 * r0 ** 2 -
                                                  L ** 3 * t0 ** 2 -
                                                  8 * L ** 2 * r0 ** 2 * t0 +
                                                  4 * L * r0 ** 2 * t0 ** 2 +
                                                  4 * r0 ** 2 * t0 ** 3)) / 2),
                         (mul, 1. / (L ** 3 - 2 * L ** 2 * t0 +
                                     L * t0 ** 2 + t0 ** 3)))
        d = (t0 - 2 * r) * L / t0 + 2. * r
        Den = np.pi * (r0 ** 2) * t0 * rho_m / (4. * r ** 2 * (t0 + L)) * 1.e12
        # Properties of formed anode
        term = alpha * Vf / P
        S = thread_first((r - term) ** 2,
                         (add, ((r - term) * (t0 - 2. * term) +
                                (d / 2. - term) *
                                (L - 2. * alpha * Vf + 2. * term))),
                         (add, ((d / 2. - term + alpha * Vf) ** 2 * (-1))),
                         (mul, 2.e-12 / (r0 * r0 * t0 * rho_m)))
        CVg = thread_first((r - term) ** 2,
                           (add, ((r - term) * (t0 - 2. * term) +
                                  (d / 2. - term) *
                                  (L - 2. * alpha * Vf + 2. * term))),
                           (add, ((d / 2. - term + alpha * Vf) ** 2 * (-1))),
                           (mul, 2. * e0 * K / (alpha * r0 * r0 * t0 * rho_m)))
        CVcc = CVg * Den
        gap = L + 2. * term * (1 - P)
        return {'r0': np.repeat(r0, len(L[d > 0])),
                't0': np.repeat(t0, len(L[d > 0])),
                'Vf': np.repeat(Vf, len(L[d > 0])),
                'L': L[d > 0],
                'r': r[d > 0],
                'd': d[d > 0],
                'Den': Den[d > 0],
                'S': S[d > 0],
                'CV/g': CVg[d > 0],
                'CV/cc': CVcc[d > 0],
                'gap': gap[d > 0]}

    def run_model(self, V):
        '''
        Run the model over a range of voltages.
        Results are stored in cap
        Parameter
            V: ndarray of voltages
        '''
        for i, Vf in enumerate(reversed(V)):
            self.cap = self._flakemodel(Vf) if i == 0 else merge_with(
                np.concatenate, self._flakemodel(Vf), self.cap)

    def get_dataV(self, V):
        '''
        After running model, get results at a given voltage.
        Parameter
            V: voltage (volts)
        Returns
            dictionary whose values for each key are ndarrays
            of model results
        '''
        return {key: self.cap[key][self.cap['Vf'] == V]
                for key in self.cap}


def compare_metals(metals_data, labels, Vf):
    '''
    Compare CV/g and CV/cc for different metals at a given formation voltage
    Parameters
    metals_data: dict with model data at V
    labels: labels for legend
    Vf: voltage corresponding to the groups voltage
    '''
    majorFormatter = FuncFormatter(
        lambda x, pos: '{:,.0f}'.format(x))
    fontsize = 18
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].set_title(''.join(['$V_f=$',  unicode(Vf), '$V$']),
                    fontsize=fontsize)
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
        e.grid(True)
    ax[0].set_ylabel('CV/g', fontsize=fontsize)
    ax[1].set_ylabel('CV/cc', fontsize=fontsize)
    plt.tight_layout()
    plt.show()
    return fig, ax


def read_props(fpath):
    '''
    Read valve metal properties from csv file.
    Parameter
        fpath: path to csv file (string)
    Returns
        vmp: dictionary of dictinaries
        {metal1: {prop1: xxx, prop2: yyy},
        metal2: {prop1: xxx, prop2: yyy}}
        where the values are floats
    '''
    with open(fpath, 'rb') as f:
        reader = csv.reader(f)
        data = [row for row in reader]
    head = data.pop(0)
    header = head[1:]
    metals = [ll.pop(0) for ll in data]
    datan = np.array([[float(datum) for datum in l] for l in data])
    properties = [dict(zip(header, l)) for l in datan]
    vmp = dict(zip(metals, properties))
    return vmp
