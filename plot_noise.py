#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys
import glob

# Publication-quality plotting
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})

def get_rms(infile):
    whdu = fits.open(infile)
    wcube = whdu[0].data

    # frequency info
    chans = np.arange(whdu[0].header["CRVAL3"], whdu[0].header["CRVAL3"]+whdu[0].header["NAXIS3"]*whdu[0].header["CDELT3"], 256*whdu[0].header["CDELT3"])
    chans = [x * 1.e-6 for x in chans]

    #undata = uncube[:,uncube.shape[1]/2,uncube.shape[2]/2]
    wdata = wcube[:,int(wcube.shape[1]/2),int(wcube.shape[2]/2)]
    nchans = wdata.shape[0]
    
    wrms = np.nanstd(wcube[:,0:50,0:50],axis=(1,2))
    
    return wrms[0]

names = ['J054829-203216', 'J082737-170020', 'J092012+21510', 'J150254-323228', 'J152146-192028', 'J161536-02554'] #sys.argv[1]
tints = [296/60, 296/60, 112/60, 296/60, 296/60, 296/60]
#ints = [41.276*60, (73*296)/60, 54*296/60, 222*296/60, 72*296/60, 37*296/60]
#ints = [10*60 + 5, (73*296)/60, 54*296/60, 10*60 + 5, 72*296/60, 37*296/60]
flux = [0.67, 1.2, 1.2, 0.67, 1.6, 3]
colors = ['red', 'blue', 'green', 'orange', 'purple', 'maroon']

print("1 sigma optical depth limits:")
for j, name in enumerate(names):
    n_obs = [2,4,8,16,32]
    files = []
    for n in n_obs:
        files.append(glob.glob(name[:7] + '*_'+str(n)+'.fits')[0])
    time = np.array(n_obs)*tints[j]

    rms = []
    for f in files:
        rms.append(get_rms(f))

    for i, n in enumerate(n_obs):
        const = np.sqrt(n) * rms[i]
        modelx = np.linspace(2, 32, 100)
        modely = const / np.sqrt(modelx)

        #plt.plot(modelx, modely, color='orange', linestyle='--')

    p = np.polyfit(np.log(time), np.log(rms), 1)
    
    n = 1

    def rms_model(t):
        return np.exp(p[1]) * t**p[0]
    time_model = np.linspace(np.min(time), ints[j])
    print(name, -np.log(1 - (rms_model(time_model)[-1] / flux[j])), '%')
    
    f = 0.1   # line/continuum ratio
    rms_desired = f * flux[j] / n
    print(rms_desired)

    #plt.axhline(rms_desired, linestyle=':', color=colors[j])
    #plt.plot(time_model, rms_model(time_model), linestyle='--', color=colors[j])
    plt.plot(time, rms, color=colors[j])
    plt.scatter(time, rms, label=name, color=colors[j])
x = np.linspace(3, 200, 100)
plt.plot(x, 1/np.sqrt(x), color='grey', linestyle='--') 
plt.ylabel('RMS noise (Jy/beam)')
plt.xlabel('Integration time (minutes)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('noise_plot_new.png', dpi=300)
plt.show()

