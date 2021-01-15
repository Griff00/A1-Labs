# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 14:22:50 2021

@author: georg
"""

from astropy.io import fits 
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo

hdulist = fits.open("./mosaic.fits")
header = hdulist[0].header

instrument_zero_point = header['MAGZPT']
instrument_zero_point_err = header['MAGZRR']

data = hdulist[0].data

#Building the histogram of data

#finding the min. and max. count values of the data
mins = []
for row in data:
    mins.append(min(row))
data_min = min(mins)

maxes = []
for row in data:
    maxes.append(max(row))
data_max = max(maxes)

#Change this to make the histogram more/less dettailed
histogram_no_of_bins = 1000

#Set up the bins
histogram_bin_size = (data_max - data_min)/histogram_no_of_bins
histogram_bin_lowest_vals = list(range(data_min,data_max,int(histogram_bin_size)))[:histogram_no_of_bins]

#Count the number of pixels in each bin
histogram_counts = np.zeros(histogram_no_of_bins)
for row in data:
    for i in row:
        bin_number = int(np.floor((i - data_min) / histogram_bin_size))
        if bin_number == histogram_no_of_bins: #The highest value causes an error
            bin_number -= 1
        histogram_counts[bin_number] += 1

#Plot- n.b the [:6] is to exclude the really high point, which is at [6]
plt.plot(histogram_bin_lowest_vals[:60],histogram_counts[:60], 'kx')
plt.xlabel('Number of counts')
plt.ylabel('Number of incidences')

#Fit Gaussian to noise

def Gaussian(x, mean, std, scaling):
    return scaling*np.e**(-(x-mean)**2/(2*std**2))

#The actual fitting to the data - settles but errors are v high- most of the vals being considered are 0. 
gaussian_params, pcov= spo.curve_fit(Gaussian,histogram_bin_lowest_vals[:60],histogram_counts[:60], p0=[3240,200,1e7])
gaussian_params_errs = np.sqrt(np.diag(pcov))

#Generating a plottable version 
g_xvals = np.linspace(data_min, data_max, 1000)
g_yvals = Gaussian(g_xvals, *gaussian_params)

plt.plot(g_xvals[:60],g_yvals[:60], 'r-')
plt.show()

print("The parameters used to calculate this fit were a mean of %f \pm %f, a standard deviation of %f \pm %f, and a scaling factor of %f \pm %f"%(gaussian_params[0], gaussian_params_errs[0], gaussian_params[1], gaussian_params_errs[1] , gaussian_params[2],gaussian_params_errs[2]))

