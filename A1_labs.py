# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 14:22:50 2021

@author: georg
"""

from astropy.io import fits 
import numpy as np
import matplotlib.pyplot as plt
hdulist = fits.open("./A1_mosaic.fits")
header = hdulist[0].header

instrument_zero_point = header['MAGZPT']
instrument_zero_point_err = header['MAGZRR']

data = hdulist[0].data

mins = []
for row in data:
    mins.append(min(row))
data_min = min(mins)

maxes = []
for row in data:
    maxes.append(max(row))
data_max = max(maxes)

histogram_no_of_bins = 100
histogram_bin_size = (data_max - data_min)/histogram_no_of_bins

histogram_bin_lowest_vals = list(range(data_min,data_max,int(histogram_bin_size)))[:histogram_no_of_bins]

histogram_counts = np.zeros(histogram_no_of_bins)

for row in data:
    for i in row:
        bin_number = int(np.floor((i - data_min) / histogram_bin_size))
        if bin_number == histogram_no_of_bins:
            bin_number -= 1
        histogram_counts[bin_number] += 1

plt.plot(histogram_bin_lowest_vals[:5] + histogram_bin_lowest_vals[6:],histogram_counts[:5] + histogram_counts[6:], 'x')
plt.xlabel('Number of counts')
plt.label('Number of incidences')
plt.show()




bin = np.zeros(100)
