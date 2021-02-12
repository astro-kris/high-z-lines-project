#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys

# Publication-quality plotting
from matplotlib import rc
#rc('text', usetex=True)
#rc('font',**{'family':'serif','serif':['serif']})

infile = sys.argv[1]
ylim = float(sys.argv[2])
#xlim = float(sys.argv[3])
whdu = fits.open(infile)
wcube = whdu[0].data

# frequency info
chans = np.arange(whdu[0].header["CRVAL3"], whdu[0].header["CRVAL3"]+whdu[0].header["NAXIS3"]*whdu[0].header["CDELT3"], 256*whdu[0].header["CDELT3"])
chans = [x * 1.e-6 for x in chans]

#undata = uncube[:,uncube.shape[1]/2,uncube.shape[2]/2]
wdata = wcube[:,int(wcube.shape[1]/2),int(wcube.shape[2]/2)]
nchans = wdata.shape[0]

wrms = np.nanstd(wcube[:,0:50,0:50],axis=(1,2))

for plotnum in range(2):
    nrows = 3
    fig = plt.figure(figsize = (20,15))
    for i in range(1,nrows+1):
        
        start =  int(nchans * (i - 1) / nrows)
        end =  int(nchans * i / nrows)
        freq_range = np.max(chans) - np.min(chans)
        freq_start = (freq_range * (i - 1) / nrows) + np.min(chans)
        freq_end = (freq_range * i / nrows) + np.min(chans)

        wdata_i = wdata[start:end]

        ax = fig.add_subplot(nrows,1,i)
        ax2 = ax.twiny()
        ax2.set_xlim((freq_start, freq_end))
        ax.fill_between(np.arange(start,end), -3.e3*wrms[start:end], 3.e3*wrms[start:end], 
                color="blue", alpha=0.3, lw=0)
        #ax.plot(undata, label="unweighted", color="red", alpha = 0.3)
        if plotnum == 0:
            ax.plot(np.arange(start,end), 1.e3*wdata_i, label="weighted", color="black")

        bands = 8
        deg = [3] * bands
        from spec_polyfit import fit_poly
        for j in range(bands):
            xaxis = np.arange(start,end)
            band_start = int(len(wdata_i) * j / bands)
            band_end = int(len(wdata_i) * (j+1) / bands)
            xraw = xaxis[band_start:band_end]
            band = wdata_i[band_start:band_end]
            x, y, poly = fit_poly(xraw, band, deg[j])

            if plotnum == 0:
                ax.plot(x, 1e3*(poly - np.mean(poly)) , label="fit", color="orange")
            else:
                ax.plot(x, 1.e3*(y - (poly - np.mean(poly))), 
                    label="weighted", color="black")

        #start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start, end, 100))
        ax2.set_xticks(np.linspace(freq_start, freq_end, 12))
        if i == 1:
            ax2.set_xlabel("Frequency / MHz")
        #ax.plot(3*wrms, ls = ":", color="black")
        #ax.plot(-3*wrms, ls = ":", color="black")
        ax.set_ylabel("peak flux density (mJy/beam)")
        if i == nrows:
            ax.set_xlabel("10-kHz channel number")
        #ax.set_ylim([-0.5, 0.5])
        ax.set_xlim([start, end])
        ax.set_ylim([-ylim, ylim])
        #f=22
        #ax.set_xlim([(128*f)+0, (128*f)+256])
        #ax.axvline(x=(128*f)+124)
        #ax.axvline(x=(128*f)+131)
        #ax.legend(loc = 2)
        #ax2.axvspan(207.06, 208.27, alpha = 0.2, color="red")
    #fig.savefig(infile.replace(".fits", "_spectrum_subtracted.pdf"), bbox_inches="tight")
    if plotnum == 0:
        ext = 'polyfit'
    else:
        ext = 'subtracted'
    fig.savefig(infile.replace(".fits", "_spectrum_" + ext + ".png"), bbox_inches="tight")

#fig.savefig("temp.png", bbox_inches="tight")

#ax.set_xlim([2500,2650])
#fig.savefig("zoom.png")
