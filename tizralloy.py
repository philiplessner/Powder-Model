# coding: utf-8
from __future__ import print_function, division, unicode_literals
import csv
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt
from powdermodel import PowderModel, compare_metals, read_props
from mplrc import mplrc

def main():
    mpl_params = mplrc()
    vmp = read_props('valve-metal-props.csv')
    
    V = np.arange(50, 250, 50)
    ad_fac = 0.38
    
    # Tantalum
    ta_anode = {'a': np.arange(0.0035, 0.2625 + 0.0035, 0.0035),
                'Ds': 6.3,
                'ratio_ptc': 0.25,
                'pore_K': 10.5,
                'lam': 0.}
    ta_model = PowderModel(ta_anode,
                           'Ta', vmp['Ta']) 
    ta_model.run_model(V)
    ta_model.plot_data()
    # TiZr
    AW_Zr = 91.2
    AW_Ti = 47.9
    AW_O = 16.0
    # Ti30Zr
    MWTi30ZrOxide = 0.7 * AW_Ti + 0.3 * AW_Zr + 2 * AW_O
    MWTi30Zr = 0.7 * AW_Ti + 0.3 * AW_Zr
    GammaTi30Zr = (4.2 * MWTi30Zr) / (5.23 * MWTi30ZrOxide)
    ti30zr_props = {'Gamma': GammaTi30Zr,
                    'Molecular Weight of metal': MWTi30Zr,
                    'Molecular Weight oxide': MWTi30ZrOxide,
                    'Moles metal/Mole oxide': 1.0,
                    'Pilling': 1. / (GammaTi30Zr * 1.0),
                    'b0': 0.0026,
                    'density of metal': 5.23e-12,
                    'density of oxide': 4.2e-12,
                    'dielectric constant': 40.,
                    'microns per volt': 0.002}
    ti30zr_anode = {'a': np.arange(0.0020, 0.262 + 0.0020, 0.0020),
                    'Ds': 2.0,
                    'ratio_ptc': 0.25,
                    'pore_K': 10.5,
                    'lam': 0.}
    ti30zr_model = PowderModel(ti30zr_anode,
                               'Ti30Zr', ti30zr_props)
    ti30zr_model.run_model(V)
    ti30zr_model.plot_data()     
    
    # Ti60Zr
    MWTi60ZrOxide = 0.4 * AW_Ti + 0.6 * AW_Zr + 2 * AW_O
    MWTi60Zr = 0.4 * AW_Ti + 0.6 * AW_Zr
    GammaTi60Zr = (4.5 * MWTi60Zr) / (5.81 * MWTi60ZrOxide) 
    ti60zr_props = {'Gamma': GammaTi60Zr,
                    'Molecular Weight of metal': MWTi60Zr,
                    'Molecular Weight oxide': MWTi60ZrOxide,
                    'Moles metal/Mole oxide': 1.0,
                    'Pilling': 1. / (GammaTi60Zr * 1.0),
                    'b0': 0.0026,
                    'density of metal': 5.81e-12,
                    'density of oxide': 4.5e-12,
                    'dielectric constant': 40.,
                    'microns per volt': 0.0019}
    ti60zr_anode = {'a': np.arange(0.0020, 0.262 + 0.0020, 0.0020),
                    'Ds': 2.21,
                    'ratio_ptc': 0.25,
                    'pore_K': 10.5,
                    'lam': 0.}
    ti60zr_model = PowderModel(ti60zr_anode,
                               'Ti60Zr', ti60zr_props)
    ti60zr_model.run_model(V)
    ti60zr_model.plot_data() 
       
    # Compare the caps at voltage
    fig, ax = compare_metals([ta_model.get_dataV(100),
                              ti30zr_model.get_dataV(100),
                              ti60zr_model.get_dataV(100)],
                             ['Ta', 'Ti30Zr', 'Ti60Zr'],
                             100)
    
    plt.clf()
    

def rolloff():
    mpl_params = mplrc()
    vmp = read_props('valve-metal-props.csv')
    
    V = np.arange(50, 160, 10)
    
    # Tantalum
    ta_anode = {'a': np.arange(0.0035, 0.2625 + 0.0035, 0.0035),
                'Ds': 6.3,
                'ratio_ptc': 0.25,
                'pore_K': 10.5,
                'lam': 0.}
    ta_model = PowderModel(ta_anode,
                           'Ta', vmp['Ta']) 
    ta_model.run_model(V)
    ta_model.make_plots()
    ta_model.v_rolloff(400.)
    
def main2():
    vmp = read_props('valve-metal-props.csv')
    
    V = np.arange(25, 175, 25)
    
    # Tantalum
    ta_anode = {'a': np.arange(0.0035, 0.2625 + 0.0035, 0.0035),
                'Ds': 6.3,
                'ratio_ptc': 0.25,
                'pore_K': 10.5,
                'lam': 0.}
    ta_model = PowderModel(ta_anode,
                           'Ta', vmp['Ta']) 
    ta_model.run_model(V)
    fpath = 'out3.csv'
    ta_model.write_csv(fpath)
    
    
if __name__ == '__main__':
    rolloff()

