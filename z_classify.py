import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
class sx():
    def __init__(self, time, distance, model, expo, mod, feature):
        self.time = time
        self.distance = distance
        self.mo = model
        self.expo = expo
        self.mod = mod
        self.feat = feature

    def av_sn(self):
        def outliers_z_score(data):
            threshold = 3
            mean = np.mean(data)
            std = np.std(data)
            z_score = [(y-mean)/std for y in data]
            return np.where(np.abs(z_score)>threshold)

        df = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.mo}/{self.mod}/{self.feat}/{self.mod}_f_{self.mo}_{self.time}_{self.distance}.csv')
        '''
        fig, ax = plt.subplots()
        ax.set_xlabel('Flux [ph cm$^{-2}$ s$^{-1}$]', fontfamily='serif')
        ax.set_ylabel('E [keV]', fontfamily='serif')
        ax.scatter(df['flux'], df['energy'])
        plt.show()
        '''
        #df = df.drop(outliers_z_score(df['energy'])[0])
        '''
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('Flux [ph cm$^{-2}$ s$^{-1}$]')
        ax.set_ylabel('Energy [keV]')
        ax.set_zlabel('S/N ')
        ax.scatter(df1['flux'], df1['energy'], df1['snr'])
        plt.show()
        '''
        
        print('Days from the explosion: ', self.time)
        print('Distance of the source: ', self.distance)
        print('FWHM: ', round(df['FWHM'].mean(),2), 'keV')
        print('Energy: ', round(df['energy'].mean(),2), 'keV')
        print('Flux: ', round(df['flux'].mean(),7), 'ph/s/cm^2')
        print('S/N: ', round(df['snr'].mean(),2))
        
        # velocity
        if df['snr'].mean() > 3:
            v = 3e5*(df['energy'].mean()-158)/df['energy'].mean()
            print('Velocity of ejecta: ', round(v,2), 'km/s')

        dt = {
                        'flux':df['flux'].mean(),
                        'energy':df['energy'].mean(),
                        'FWHM':df['FWHM'].mean(),
                        'snr':df['snr'].mean(),
                        'distance':self.distance,
                }

        df_d = pd.DataFrame(data=dt, index=[1])
        #df_d.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}_norm/{self.mo}/{self.mod}/{self.feat}/{self.mod}_pippo_{self.mo}_{self.time}.csv',index=False)

        try:
            df1_d = pd.read_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.mo}/{self.mod}/{self.feat}/{self.mod}_pippo_{self.mo}_{self.time}.csv')
            df2_d = pd.concat([df_d, df1_d])
            df2_d.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.mo}/{self.mod}/{self.feat}/{self.mod}_pippo_{self.mo}_{self.time}.csv',index=False)
        except FileNotFoundError:
            df_d.to_csv(f'/mnt/c/_Home/Type_Ia_SN_simulations_ASTENA/expo_{self.expo}/{self.mo}/{self.mod}/{self.feat}/{self.mod}_pippo_{self.mo}_{self.time}.csv',index=False)

