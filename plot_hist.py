#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from astropy.io import fits
import sys
import glob

# Publication-quality plotting
from matplotlib import rc
#rc('text', usetex=True)
#rc('font',**{'family':'serif','serif':['serif']})

def plot_hist(filenames):
	#filenames = sys.argv[1]
	#print(filenames)
	infiles = glob.glob(filenames)
	for i, f in enumerate(infiles):
		whdu = fits.open(f)
		wcube = whdu[0].data

		wdata = wcube[:,int(wcube.shape[1]/2),int(wcube.shape[2]/2)]
		nchans = wdata.shape[0]

		wrms = np.nanstd(wcube[:,0:50,0:50],axis=(1,2))

		wdata = wdata[wdata != 0]
		std = np.std(wdata)
		mean = np.mean(wdata)
		h = 3.5 * std / nchans**(1/3)
		nbins = int((np.max(wdata) - np.min(wdata)) / h)

		x = np.linspace(np.min(wdata), np.max(wdata), 100)
		
		plt.subplot(len(infiles),1,i+1)
		plt.hist(wdata, nbins, density=True)
		plt.plot(x, norm.pdf(x, mean, std))
	plt.show()
	
