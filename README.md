# A1-Labs
A1 Labs Image Processing Project

This branch contains 'pixel-by-pixel' magntiude method. See double aperture branch for 'double aperture method' (recommended).

This branch contains 2 scripts:

Detect_whole.py - This script contains all the functionality related to Reading the FITS data, plotting and calculating the photon count histogram and gaussian fit (can be toggled on/off using boolean variable). Source detection algorithm, mask and contamination removal functionalities are included. Magnitude calculations are performed using 'pixel-by-pixel' method. Output data is saved as cdv. Requires presence of FITS file.

N_vs_m_plot.py - This script reads in galaxy_list.csv and plots the log(N) vs m plot. Error propagation and fitting functionality is included here. Requires presence of CSV file.

Standard packages are required: numpy, scipy, matplotlib, astropy
