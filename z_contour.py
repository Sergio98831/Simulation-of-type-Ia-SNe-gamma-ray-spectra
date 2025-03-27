from heasp import *
from xspec import *
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib import cm
from matplotlib import colors
from scipy import optimize
from matplotlib.colors import LogNorm
font = {'family': 'serif',
        'weight': 'normal',
        }

class CONTOUR():
    def __init__(self, model, distance, time, expo, mode, feature):
        self.model = model
        self.distance = distance
        self.time = time
        self.expo = expo
        self.mod = mode
        self.feat = feature

    def pippo(self):
        snr, s, d, s1 = [], [], [], []
        for t in self.time:
            df = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.model}/{self.mod}/{self.feat}/{self.mod}_ugo_{self.model}_{round(t,2)}.csv')
            df1 = df.iloc[::-1]
            s.append(df1['snr'])
        snr = np.array(s)
        fig, ax = plt.subplots()
        levels1 = np.linspace(0,10,11)
        levels2 = [3,5]
        CS = ax.contourf(df1['distance'], self.time, snr, levels1)
        cs = ax.contour(df1['distance'], self.time, snr, levels2, colors=['y','r'])
        cb = fig.colorbar(CS)
        #ax.hlines(y = 17.5, xmin = min(self.distance), xmax = 3.5, linewidth=0.5, linestyles='dashed', color='w')
        #ax.vlines(x = 3.5, ymin = min(self.time), ymax = 17.5, linewidth=0.5, linestyles='dashed', color='w')
        ax.scatter(3.5,17.5, c= 'orange', label='SN2014J')
        ax.scatter(6.4,20.0, c='w', label='SN2011fe')
        #ax.text(3.0, 18.0, 'SN2014J', color='k', fontdict=font)
        cb.set_label('S/N Ratio', fontdict = font)
        ax.set_xlabel('Distance [Mpc]', fontdict=font, fontsize=10)
        ax.set_ylabel('Time from explosion [days]', fontdict=font, fontsize=10)
        ax.set_xlim(min(df['distance']), max(df['distance']))
        ax.legend()
        plt.savefig(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/ASTENA/{self.mod}_contour_sigma_{self.model}_{self.expo}.pdf')
        plt.show()
        plt.close()
        


    def sn(self):

        def func(x, b, c):
            y = b/(x) + c
            return y

        for t in self.time:
            df = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.model}/{self.mod}/{self.feat}/{self.mod}_pippo_{self.model}_{round(t,2)}.csv', header=0)
            x, y, fortitudo, energy, FW = [], [], [], [], []
            for element, sigma, e, flux in zip(df['distance'], df['snr'], df['energy'], df['flux']):
                x.append(element)
                fortitudo.append(1/element)
                y.append(sigma)
                energy.append(e)
                #FW.append(fwhm)

            a,b = np.polyfit(fortitudo, y, 1)
            #fig, (ax, ax1) = plt.subplots(2)
            gif, ax = plt.subplots()
            xx = np.linspace(min(x), max(x), 1000)
            yy = func(xx, a,b)
            x1 = np.asarray(x)
            y1 = func(x1, a, b)
            res = y - y1
            sq = res**2
            standard_error = (sum(sq)/(len(y)-2))**0.5
            #s_d_p = standard_error/b
            ax.set_ylabel('S/N', font)
            ax.set_xlabel('Distance [Mpc]', font)
            ax.loglog(x, y, 'k.')
            ax.set_ylim(0.5,max(y)+5)
            ax.set_xlim(0.5,max(self.distance) + 5)
            #ax.grid(True)
            popt, pcov = optimize.curve_fit(func, xx, yy)
            ax.loglog(xx, func(xx, *popt), color = 'red')
            f = open(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/ASTENA/prova_{self.model}_{t}_{self.expo}.txt', 'w')
            f.write('------------------------------------------------------\n')
            f.write(f'Days from the explosion: {t}\n')
            f.write(f'Centroid: {round(np.mean(e), 1)} keV\n')
            #f.write(f'Flux: {np.mean(flux)} ph/cm^2/s\n')
            #f.write(f'FWHM: {np.mean(FW)} keV')
            f.write(f'Velocity: {round(3e5*(np.mean(energy)-158)/np.mean(energy), 1)}\n')
            f.write(f'At 3.5 Mpc: {round(func(3.5, *popt),1)} \u03C3 \n')
            f.write(f'At 6.4 Mpc: {round(func(6.4, *popt),1)} \u03C3 \n')
            num = [5,3]
            col = ['b', 'g']
            for i,j in zip(num,col):
                puma = np.where(yy>i)
                c = np.asarray(puma)
                if c.size !=0:
                    f.write(f'{i}\u03C3 at {round(xx[c[0][-1:]][0],1)} Mpc\n')

                    ax.loglog(xx[c[0][-1:]], yy[c[0][-1:]], '^', label=f'{i}\u03C3 at {round(xx[c[0][-1:]][0],1)} Mpc', color=j)
                    ax.vlines(x = xx[c[0][-1:]], ymin=0, ymax=yy[c[0][-1:]], color='gray', linestyle=':', linewidth=0.5)
                    ax.hlines(y = yy[c[0][-1:]], xmin=0, xmax=xx[c[0][-1:]], color='gray', linestyle=':', linewidth=0.5)
                else:
                    pass
            f.close()
            ax.legend()

            plt.savefig(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/ASTENA/{self.mod}_Sigma_{self.model}_at_{round(t,2)}_days_{self.expo}.pdf')
            plt.close()
            ds = {
                'distance':xx,
                'snr':yy
            }
            df1 = pd.DataFrame(ds)
            df1.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.model}/{self.mod}/{self.feat}/{self.mod}_ugo_{self.model}_{round(t,2)}.csv',index=False)

            os.remove(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.model}/{self.mod}/{self.feat}/{self.mod}_pippo_{self.model}_{round(t,2)}.csv')

class xcount():
    def __init__(self, model, distance, time, expo, mode):
        self.model = model
        self.distance = distance
        self.time = time
        self.expo = expo
        self.mod = mode

    def p(self):

        def choise(p1, p2, ax1, ax2):
            for ax, p in zip((ax1, ax2), (p1, p2)):
                levels = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
                CS = ax.contourf(self.distance, (self.time), p, levels, norm = LogNorm())
                cs = ax.contour(self.distance, self.time, p, [0.01], colors=['r'], norm = LogNorm())
                ax.vlines(x = 3.5, ymin = min(self.time), ymax = max(self.time), linewidth=0.5, linestyles='dashed', color='k')
                ax.vlines(x = 6.4, ymin = min(self.time), ymax = max(self.time), linewidth=0.5, linestyles='dashed', color='k')
                ax.scatter(3.5, 17.5,color='orange', label='SN2014J')
                ax.scatter(6.4, 20.0, color='k', label='SN2011fe')
                ax.set_xlabel('Distance [Mpc]', fontdict = font)
                cb = fig.colorbar(CS, location='top', orientation='horizontal')
                cb.set_label('P-Value', fontdict = font)
                ax.legend()
            plt.savefig(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/ASTENA/{self.mod}_contour_p_value_{self.model}_{self.expo}.pdf')
            plt.show()
            plt.close()

        p_ddt_tot, p_d2t_tot, p_merge_tot = [], [], []
        for z in self.time:
            df = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.model}/{self.mod}_pval_{self.model}_{z}_{self.expo}.csv')
            p_ddt_tot.append(df['ddt'])
            p_d2t_tot.append(df['d2t'])
            p_merge_tot.append(df['merge'])

        fig, (ax1, ax2) = plt.subplots(ncols=2)
        ax1.set_ylabel('Time [days]', fontdict=font)

        #fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap='magma'), ax=ax2, orientation='vertical', label='a colorbar label')
        if self.model == 'ddt':
            choise(p_d2t_tot, p_merge_tot, ax1, ax2)
        elif self.model == 'd2t':
            choise(p_ddt_tot, p_merge_tot, ax1, ax2)
        else:
            choise(p_ddt_tot, p_d2t_tot, ax1, ax2)
