
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
hdulist = fits.open("mosaic.fits")

#Test data
data = [[10,1,1,2,1,2,1,2,1,2,1,10],
        [1,10,1,2,1,10,1,2,1,2,1,10],
        [1,1,1,2,10,10,10,2,1,2,1,2],
        [10,2,1,2,1,2,1,2,1,2,1,10]]
dimensions = (len(data),len(data[0]))

#Mask is altered when sources are found. Sources are given integer identifiers which are added into mask.
mask = np.zeros(dimensions)
mean = np.mean(data)
background = mean #Just tempory!! Until we have actual background

def source_detection(data,mask,background,dimensions): 
    # Cyles through each pixel checking for sources
    count = 0
    for y in range(0,len(data)):
        for x in range(0,len(data[0])):
            if data[y][x] > background:
                count, mask = source_handler(data,x,y,mask,background,dimensions,count)
    return count,mask
                
def source_handler(data,x,y,mask,background,dimensions,count):
    # Once a source is found, need to determine if pixel is part of a mulit-pixel distributed source
    x_low = x-1
    x_high = x+1
    y_high = y+1
    y_low = y-1

    if x_low < 0:
        x_low = 0
    if x_high > dimensions[1]:
        x_high = dimensions[1]
    if y_low < 0:
        y_low = 0
    if y_high > dimensions[0]:
        y_high = dimensions[0]                    

    for j in range(y_low,y_high):         
        for i in range(x_low,x_high):
           # print("here",x,y)
            #print("i",i,"j",j,"mask",mask[j][i])
            if mask[j][i] > 0:
                mask[y][x] = mask[j][i]         
    if  mask[y][x] == 0:
        count +=1
        mask[y][x] = count
    return count, mask

count, mask = source_detection(data,mask,background,dimensions)
print(mask,count)