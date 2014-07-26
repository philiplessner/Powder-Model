# coding: utf-8
from __future__ import print_function, division, unicode_literals
from itertools import chain
from IPython.display import HTML
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from powdermodel import PowderModel
from utils import format_HTML


class PowderModelPandas(PowderModel):

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
