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

print(data)
