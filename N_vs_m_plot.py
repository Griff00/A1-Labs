# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 12:29:40 2021

@author: Jack
"""


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('galaxy_list.csv', delimiter=",")

count = len(data)
data_min = min(data)

data_max = max(data)
histogram_no_of_bins = 50
histogram_bin_size = (data_max - data_min)/histogram_no_of_bins
histogram_bin_lowest_vals = np.linspace(data_min,data_max,histogram_no_of_bins)
histogram_counts = np.zeros(histogram_no_of_bins)

for i in data:
    bin_number = int(np.floor((i - data_min) / histogram_bin_size))
    if bin_number == histogram_no_of_bins: #The highest value causes an error
        bin_number -= 1
    histogram_counts[bin_number] += 1

N_cum = np.cumsum(histogram_counts)

plt.plot(histogram_bin_lowest_vals,N_cum)
plt.plot(histogram_bin_lowest_vals,N_cum,'x')
plt.xlabel('Magnitude (m)')
plt.ylabel('N(<m)')
plt.show()

log_N = np.log10(N_cum)

plt.figure()
plt.plot(histogram_bin_lowest_vals,log_N)
plt.plot(histogram_bin_lowest_vals,log_N,'x')
plt.xlabel('Magnitude (m)')
plt.ylabel('log(N(<m)')
plt.show()

print(count)