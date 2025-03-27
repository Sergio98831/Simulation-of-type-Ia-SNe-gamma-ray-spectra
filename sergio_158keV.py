import numpy as np
from xspec import *
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import os
import pandas as pd
import time as tt
from z_classify import sx
from z_contour import CONTOUR
snr, snrr, snrrr, ff = [], [], [], []
distance = np.arange(1.0,11.0,1.0)
model = ['ddt', 'd2t', 'merge']
expo = float(input('Please insert the Exposure Time: '))
mode = input('"old" or "new" rmf: ')
old = {
    'rmf':'astena-nft.rmf',
    'arf':'astena-nft.arf',
    'bkg':'astena-nft.bkg'
}

new = {
    'rmf':'new-astena-nft.rmf',
    'arf':'new-astena-nft.arf',
    'bkg':'astena-nft.bkg'
}

font = {
        'family': 'serif',
        'weight': 'normal',
        'size': 15,
        }

dt2 = expo/(24*3600) # converts exposure time in days
time = np.arange(10.0, 61.0, 5.0)

c, d, FWHM, g, aq, aq1 = [], [], [], [], [], [] # features of emission lines
f = {
        'ddt':0.45,
        'd2t':0.30,
        'merge':0.15
        }

num = int(input('How many sources do you want to simulate? '))
feature = '158keV'
start = tt.process_time()
for milf in model:
        for x in time:
                for i in distance:
                        for l in range(num):
                                norm = 1/(i**2)
                                x_time = (x + x + dt2)/2
                                Xset.chatter = 5
                                AllData.dummyrsp(50,4000) # dummy response
                                # simulations: selecting the model and freezing the time parameter
                                if milf == 'ddt':
                                        m = Model('atable{ddt.mod}',setPars={1:x_time, 2:norm})
                                        m.ddt.time.frozen = True
                                elif milf == 'd2t':
                                        m = Model('atable{d2t.mod}',setPars={1:x_time, 2:norm})
                                        m.d2t.time.frozen = True
                                else:
                                        m = Model('atable{merge.mod}', setPars={1:x_time, 2:norm})
                                        m.merge.time.frozen = True
                                # selecting the response matrix
                                if mode == 'old':
                                        fs1 = FakeitSettings(response=old['rmf'], arf=old['arf'], background=old['bkg'], exposure=expo, correction=1, backExposure=expo, fileName=f'sn_{l}.fak')
                                        resp = old['rmf']
                                elif mode == 'new':
                                        fs1 = FakeitSettings(response=new['rmf'], arf=new['arf'], background=new['bkg'], exposure=expo, correction=1, backExposure=expo, fileName=f'sn_{l}.fak')
                                        resp = new['rmf']
                                else:
                                        ValueError

                                # simulations
                                fake_file = AllData.fakeit(1, fs1)
                                os.system(f'ftgrouppha infile=sn_{l}.fak backfile=sn_{l}_bkg.fak outfile=sn_{l}_grp25.pha grouptype=min groupscale=50 respfile={resp}')
                                AllData.clear()

                                # evaluation of S/N for 158 keV feature
                                
                                file = f'sn_{l}_grp25.pha'
                                s = Spectrum(file) # spectrum of the source      
                                Plot.xAxis = "keV" # configures the x axis as "energy"
                                s.ignore("**-120.0 250.0-**") # limit our attention around the emission line energy
                                Plot.device = "/null" # I imported matplotlib
                                Fit.perform() # performs the chi2 fit
                                Plot("ufspec") # save flux and energy
                                xVals1 = np.array(Plot.x()) # x = energy
                                yVals1 = np.array(Plot.y()) # y = flux
                                v = np.array(Plot.model()) # value of the model (theoretical flux-y value)
                                x1_err = np.array(Plot.xErr()) # evaluation of energy and flux errors
                                y1_err = np.array(Plot.yErr())

                                # select the counts to evaluate S/N

                                Plot('counts')
                                x1 = np.array(Plot.x()) # energy
                                y1 = np.array(Plot.y()) # simulated counts
                                xErr1 = np.array(Plot.xErr())
                                yErr1 = np.array(Plot.yErr())
                                AllModels.calcFlux('155 170 err') # flux within the 155-170 keV region
                                flux = s.flux
                                c.append(flux[3])
                                e = abs((flux[4]-flux[5])/2)

                                # spectrum of the continuum
                                s.ignore(f"1:160.0-170.0") # ignore the emission line      
                                AllModels.clear() # removing the simulation model
                                m4 = Model("pegpw", setPars={1:1.0, 2:120, 3:250, 4:flux[0]}) # the spectra can be describe with a power law
                                Fit.query = 'yes'
                                Fit.perform() # we find the continuum model
                                a = m4.pegpwrlw.PhoIndex.values[0]
                                b = m4.pegpwrlw.norm.values[0]
                                s.notice('all')
                                s.ignore("**-120.0 250.0-**")
                                Plot('ufspec') # we plot the continuum model considering the emission line range
                                xVals3 = np.array(Plot.x())
                                yVals3 = np.array(Plot.model())
                                Plot("counts")
                                x3 = np.array(Plot.x())
                                y3 = np.array(Plot.model())
                                yErr3 = np.sqrt(y3)
                                AllModels.clear()

                                #Gauss -> to find the centroid and FWHM
                                s.notice('all')
                                s.ignore("**-100.0 250.0-**")
                                m5 = Model("pegpw + gauss")
                                m5.setPars(a,120, 250, b, 164, 3, 1e-5)
                                m5.pegpwrlw.PhoIndex.frozen = True # freezing powerlaw parameter
                                m5.pegpwrlw.norm.frozen = True
                                Fit.perform()
                                d.append(m5.gaussian.LineE.values[0]) # centroid of the line
                                aq.append(m5.gaussian.LineE.sigma) # error on the centroid
                                aq1.append(m5.gaussian.Sigma.sigma) # FWHM error
                                FWHM.append(m5.gaussian.Sigma.values[0]) # FWHM
                                Plot("ufspec")
                                xVals2 = np.array(Plot.x())
                                yVals2 = np.array(Plot.model())
                                yErr2 = np.sqrt(yVals2)

                                # to evaluate the S/N we need to integrate a function -> interpolate the counts to obtain an integrable function
                                n = 100
                                xfine = np.linspace(155, 170, n*10) # integration range

                                y0 = interpolate.interp1d(x1, y1, kind='quadratic', fill_value='extrapolate') # emission line
                                z0 = interpolate.interp1d(x3, y3, kind='quadratic', fill_value='extrapolate') # continuum

                                # integrate
                                Iy0 = integrate.simpson(y0(xfine), x=xfine) # total counts
                                Iz0 = integrate.simpson(z0(xfine), x=xfine) # continuum counts

                                if Iy0 > 0 and Iz0 > 0: # if they are > 0 we can evaluate:
                                        Cl = Iy0 - Iz0  # emission line counts
                                        I_err = np.sqrt(Iy0) 
                                        Z_err = np.sqrt(Iz0)
                                        noise = np.sqrt((I_err)**2 + (Z_err)**2) # total errors
                                        q = Cl/Iy0
                                        snr.append(Cl/noise)
                                else:
                                        snr.append(0)
                                        q = 0

                                g.append(q)
                                snrr = np.array(snr)
                                snrrr = np.nan_to_num(snrr)

                                ff.append(np.mean(snrrr))

                                '''fig, ax = plt.subplots()
                                ax.set_xlabel('Energy [keV]', fontfamily='serif')
                                ax.set_ylabel('Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$', fontfamily='serif')
                                ax.set_xlim(120, 220)
                                ax.errorbar(xVals1, yVals1, abs(y1_err), x1_err, ',',linewidth=0.4, color='black')

                                #ax.plot(xVals1, v, linewidth=1.5, color='r')
                                ax.plot(xVals2, yVals2, linewidth=1.5, color='r')
                                ax.plot(xVals3, yVals3, linewidth=1.5, linestyle='-.', color ='b')
                                plt.savefig(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/SIMULATIONS_teo/{mode}_z_{milf}_{l}_{x}_days_{i}_Mpc.jpeg')
                                plt.show()
                                plt.close()'''

                                AllData.clear()
                                AllModels.clear()
                                os.system('rm sn_*')
                                snr = []

                        dt = {
                                'flux':c,
                                'energy':d,
                                'sigma_e':aq,
                                'FWHM':FWHM,
                                'sigma_f':aq1,
                                'snr':ff,
                                'f':g,
                                'distance':np.full_like(ff,i),
                                'time':np.full_like(ff,x)
                        }

                        df = pd.DataFrame(data=dt)

                        try:
                                df1 = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/SIMULATIONS_teo/{expo}_{mode}_f_{milf}_{x}_{i}.csv')
                                df2 = pd.concat([df, df1])
                                df2.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/SIMULATIONS_teo/{expo}_{mode}_f_{milf}_{x}_{i}.csv',index=False)
                        except FileNotFoundError:
                                df.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/SIMULATIONS_teo/{expo}_{mode}_f_{milf}_{x}_{i}.csv',index=False)

                        sx(round(x,2), i, milf, expo, mode, feature).av_sn() # preparing data

                        c, d, ff, FWHM, g, aq, aq1 = [], [], [], [], [], [], []

        CONTOUR(milf, distance, time, expo, mode, feature).sn() # significance of 158-keV line as a function of distance
        CONTOUR(milf, distance, time, expo, mode, feature).pippo() # significance as a function of distance and time

end = tt.process_time()
secondi = int(end - start)
ss      = secondi %  60
minuti  = secondi // 60
mm      = minuti  %  60
hh      = minuti  // 60

print(secondi, "s ->", hh, "h", mm, "m", ss, "s")
