import numpy as np
from heasp import *
from xspec import *
import matplotlib.pyplot as plt
from scipy import interpolate, integrate, signal
import os
import pandas as pd
import time as tt
from z_classify import sx
from z_contour import CONTOUR
snr, snrr, snrrr, ff = [], [], [], []
distance = np.arange(1.0,22.0,2.0)
#distance = [3.5]
#model = ['ddt', 'd2t', 'merge']
model = ['ddt']
expo = float(input('Please insert the Exposure Time: '))
mode = input('Do you want old or new response matrix? Type "old" or "new": ')
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

dt2 = expo/(24*3600)
#time = np.arange(15.0, 30.0, dt)
#time = np.arange(10.0, 61.0, 5.0)
time = [17.5]
font = {
        'family': 'serif',
        'weight': 'normal',
        }
c, d, FWHM, g, aq, aq1 = [], [], [], [], [], []
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
                                AllData.dummyrsp(50,4000) #dummy response
                                #simulations
                                if milf == 'ddt':
                                        m = Model('atable{ddt.mod}',setPars={1:x_time, 2:norm})
                                        m.ddt.time.frozen = True
                                elif milf == 'd2t':
                                        m = Model('atable{d2t.mod}',setPars={1:x_time, 2:norm})
                                        m.d2t.time.frozen = True
                                else:
                                        m = Model('atable{merge.mod}', setPars={1:x_time, 2:norm})
                                        m.merge.time.frozen = True

                                if mode == 'old':
                                        fs1 = FakeitSettings(response=old['rmf'], arf=old['arf'], background=old['bkg'], exposure=expo, correction=1, backExposure=expo, fileName=f'sn_{l}.fak')
                                        resp = old['rmf']
                                elif mode == 'new':
                                        fs1 = FakeitSettings(response=new['rmf'], arf=new['arf'], background=new['bkg'], exposure=expo, correction=1, backExposure=expo, fileName=f'sn_{l}.fak')
                                        resp = new['rmf']
                                else:
                                        ValueError

                                fake_file = AllData.fakeit(1, fs1)
                                os.system(f'ftgrouppha infile=sn_{l}.fak backfile=sn_{l}_bkg.fak outfile=sn_{l}_grp25.pha grouptype=min groupscale=50 respfile={resp}')
                                AllData.clear()

                                # evaluation of S/N for 158 keV feature
                                file = f'sn_{l}_grp25.pha'
                                #file = f'sn_{l}.fak'
                                s = Spectrum(file)      # spectrum of the source
                                Plot.xAxis = "keV"
                                s.ignore("**-120.0 250.0-**")
                                Plot.device = "/null"
                                Fit.perform()
                                Plot("ufspec")
                                xVals1 = np.array(Plot.x())
                                yVals1 = np.array(Plot.y())
                                v = np.array(Plot.model())
                                x1_err = np.array(Plot.xErr())
                                y1_err = np.array(Plot.yErr())
                                p = np.argmax(v)
                                peaks = [p]
                                uv = xVals1[p]
                                Plot('counts')
                                x1 = np.array(Plot.x())
                                y1 = np.array(Plot.y())
                                xErr1 = np.array(Plot.xErr())
                                yErr1 = np.array(Plot.yErr())
                                AllModels.calcFlux('155 170 err')
                                flux = s.flux
                                c.append(flux[3])
                                #d.append(uv)
                                e = abs((flux[4]-flux[5])/2)

                                s.ignore(f"1:160.0-170.0")      # spectrum of the continuum
                                Plot("ufspec")
                                xVals4 = np.array(Plot.x())
                                yVals4 = np.array(Plot.model())
                                x4_err = np.array(Plot.xErr())
                                y4_err = np.array(Plot.yErr())
                                Plot('counts')
                                x4 = np.array(Plot.x())
                                y4 = np.array(Plot.y())
                                xErr4 = np.array(Plot.xErr())
                                yErr4 = np.array(Plot.yErr())

                                AllModels.clear()
                                #s.notice('all')
                                #s.ignore("**-120.0 250.0-**")
                                m4 = Model("pegpw", setPars={1:1.0, 2:120, 3:250, 4:flux[0]})
                                Fit.query = 'yes'
                                Fit.perform()
                                a = m4.pegpwrlw.PhoIndex.values[0]
                                b = m4.pegpwrlw.norm.values[0]
                                s.notice('all')
                                s.ignore("**-120.0 250.0-**")
                                Plot('ufspec')
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
                                m5.pegpwrlw.PhoIndex.frozen = True
                                m5.pegpwrlw.norm.frozen = True
                                #m5.gaussian.Sigma.frozen = True
                                Fit.perform()
                                #print(m5.gaussian.norm.values[0])
                                #print(m5.pegpwrlw.norm.values[0])
                                #print(m5.gaussian.norm.values[0]+ m5.pegpwrlw.norm.values[0])
                                d.append(m5.gaussian.LineE.values[0])
                                aq.append(m5.gaussian.LineE.sigma)
                                aq1.append(m5.gaussian.Sigma.sigma)
                                FWHM.append(m5.gaussian.Sigma.values[0])
                                Plot("ufspec")
                                xVals2 = np.array(Plot.x())
                                yVals2 = np.array(Plot.model())
                                yErr2 = np.sqrt(yVals2)

                                n = 100
                                #xfine = np.linspace(m5.gaussian.LineE.values[0]-m5.gaussian.Sigma.values[0], m5.gaussian.LineE.values[0]+m5.gaussian.Sigma.values[0], n*10)
                                xfine = np.linspace(155, 170, n*10)
                                cfine = np.linspace(120, 250, n*10)

                                y0 = interpolate.interp1d(x1, y1, kind='quadratic', fill_value='extrapolate')
                                z0 = interpolate.interp1d(x3, y3, kind='quadratic', fill_value='extrapolate')
                                #x0 = interpolate.interp1d(xVals2, yVals2, kind='quadratic', fill_value='extrapolate')

                                #integrate
                                Iy0 = integrate.simpson(y0(xfine), x=xfine)
                                Iz0 = integrate.simpson(z0(xfine), x=xfine)

                                Cl = Iy0 - Iz0
                                if Iy0 > 0 and Iz0 > 0:
                                        Cl = Iy0 - Iz0
                                        I_err = np.sqrt(Iy0)
                                        Z_err = np.sqrt(Iz0)
                                        noise = np.sqrt((I_err)**2 + (Z_err)**2)
                                        if Cl > 0:
                                                q = Cl/Iy0
                                                snr.append(Cl/noise)
                                        else:
                                                q = 0
                                                snr.append(0)

                                else:
                                        snr.append(0)
                                        q = 0

                                g.append(q)
                                snrr = np.array(snr)
                                snrrr = np.nan_to_num(snrr)

                                ff.append(np.mean(snrrr))

                                '''fig, ax = plt.subplots()
                                #ax.set_title(rf'{milf} at {i} Mpc')
                                ax.set_xlabel('Energy [keV]', fontfamily='serif')
                                ax.set_ylabel('Photons cm$^{-2}$ s$^{-1}$ keV$^{-1}$', fontfamily='serif')
                                ax.set_xlim(120, 220)
                                ax.errorbar(xVals1, yVals1, abs(y1_err), x1_err, ',',linewidth=0.4, color='black')

                                ax.plot(xVals1, v, linewidth=1.5, color='r')
                                ax.plot(xVals2, yVals2, linewidth=1.5, linestyle='--', color='b')
                                ax.plot(xVals3, yVals3, linewidth=1.5, linestyle='-.', color ='g')
                                #ax.plot(cfine, x0(cfine), linewidth=0.75, color='k')
                                #ax.vlines(uv - 6, -4e-5, 4e-5, colors='k', linestyles='-.', linewidth=0.7)
                                #ax.vlines(uv + 6, -4e-5, 4e-5, colors='k', linestyles='-.', linewidth=0.7)
                                #ax.vlines(155, -5e-6, 5e-6, colors='k', linestyles='-.', linewidth=0.7)
                                #ax.vlines(170, -5e-6, 5e-6, colors='k', linestyles='-.', linewidth=0.7)
                                #ax.plot(xVals1[peaks], v[peaks], 'x', color='y')

                                plt.savefig(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{expo}/{milf}/{mode}_z_{milf}_{l}_{x}_days_{i}_Mpc.pdf')
                                plt.show()
                                #plt.close()'''

                                '''sn_project = {
                                        'energy': xVals1,
                                        'flux': yVals1
                                }
                                ddf = pd.DataFrame(data=sn_project)
                                ddf.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{expo}/{milf}/spectra/{mode}_sp_{l}_{milf}_{round(x,2)}_{i}.csv')
                                '''
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
                        #df.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{expo}/{milf}/{mode}/{feature}/{mode}_f_{milf}_{round(x,2)}_{i}.csv',index=False)

                        try:
                                df1 = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{expo}/{milf}/{mode}/{feature}/{mode}_f_{milf}_{x}_{i}.csv')
                                df2 = pd.concat([df, df1])
                                df2.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{expo}/{milf}/{mode}/{feature}/{mode}_f_{milf}_{x}_{i}.csv',index=False)
                        except FileNotFoundError:
                                df.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{expo}/{milf}/{mode}/{feature}/{mode}_f_{milf}_{x}_{i}.csv',index=False)

                        sx(round(x,2), i, milf, expo, mode, feature).av_sn()

                        c, d, ff, FWHM, g, aq, aq1 = [], [], [], [], [], [], []

        CONTOUR(milf, distance, time, expo, mode, feature).sn()
        #CONTOUR(milf, distance, time, expo, mode, feature).pippo()

end = tt.process_time()
secondi = int(end - start)
ss      = secondi %  60
minuti  = secondi // 60
mm      = minuti  %  60
hh      = minuti  // 60

print(secondi, "s ->", hh, "h", mm, "m", ss, "s")
