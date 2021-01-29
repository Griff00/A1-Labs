# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 12:29:40 2021

@author: Jack
"""


from astropy.io import fits
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt

input_data = np.loadtxt('galaxy_list.csv', delimiter=",")
data = input_data[0]
data_errs = input_data[1]

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

cum_counts_errs = []
for i in range (0,len(N_cum)):
    cum_counts_errs.append(np.sqrt(N_cum[i]))
    

plt.plot(histogram_bin_lowest_vals,N_cum)
plt.errorbar(histogram_bin_lowest_vals,N_cum,yerr = cum_counts_errs,fmt ='.',capsize=3)
plt.xlabel('Magnitude (m)')
plt.ylabel('N(<m)')
plt.show()

log_N = np.log10(N_cum)

log_errs = []
for i in range (0,len(log_N)):
    err = 0.434*(cum_counts_errs[i]/N_cum[i])
    log_errs.append(err)


def theory(x,grad,const):
    return grad*x + const
 #Can change values straight line is fit too (atm its 15-38)
theory_params, pcov= spo.curve_fit(theory,histogram_bin_lowest_vals[15:38],log_N[15:38], p0=[0.6,0])
theory_params_errs = np.sqrt(np.diag(pcov))

t_xvals = np.linspace(data_min, data_max, 1000)
t_yvals = theory(t_xvals, *theory_params)

plt.figure()
#plt.plot(histogram_bin_lowest_vals,log_N)
plt.errorbar(histogram_bin_lowest_vals,log_N,yerr = log_errs,fmt ='x',capsize=3)
plt.plot(t_xvals,t_yvals,label = "Straight line fit")
plt.xlabel('Magnitude (m)')
plt.ylabel('log(N(<m)')
plt.show()



graphical_params = {
  # 'axes.titlesize': 'x-large', 
  # 'axes.labelsize': 'large',sss
   'font.size': 14,
  # 'legend.fontsize': 14,
  # 'xtick.labelsize': 10,
  # 'ytick.labelsize': 10,
   'figure.figsize': [11.69,8.27],
   'errorbar.capsize': 8
   } 
plt.rcParams.update(graphical_params)
plt.grid()
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth='0.5', alpha = 0.7, color='black')
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='grey')
#plt.xlim(0,4000)
#plt.ylim(0,30000)
#plt.xticks(np.arange(1,L+1,1))
#plt.title("Height of pile with time for Oslo Model")
#plt.xlabel('Time (in number of grains added)')
#plt.ylabel('Height of pile') 
plt.style.use('default')
plt.legend;
plt.show()



print(count)
print("The parameters used to calculate this fit were a gradient of %f \pm %f, and a constant of %f \pm %f"%(theory_params[0], theory_params_errs[0], theory_params[1], theory_params_errs[1]))