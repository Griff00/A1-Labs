
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
hdulist = fits.open("mosaic.fits")

#Test data
data = [[1,1,1,2,1,2,1,2,1,2,1,1],
        [1,1,1,2,1,30,30,2,1,2,40,1],
        [1,1,1,2,1,30,30,2,1,2,1,1],
        [1,2,1,2,1,2,1,2,1,2,1,1]]
dimensions = (len(data),len(data[0]))

#Mask is altered when sources are found. Sources are given integer identifiers which are added into mask.

mean = np.mean(data)
std = np.std(data)
background = mean + std #Just tempory!! Until we have actual background

mask = np.ones(dimensions, dtype=bool)
source_map = np.zeros(dimensions)

frame_x_left = 1
frame_x_right = 1
frame_y_top = 1
frame_y_bottom = 1
frame_cutoffs = (frame_x_left,frame_x_right,frame_y_top,frame_y_bottom)


class Source:

    source_count = 0
    source_list = []
   
    def __init__(self,identifier,location,count):
                 
        self.identifier = identifier
        self.pixel_coords = []
        self.pixel_counts = []
        self.pixel_coords.append(location)
        self.pixel_counts.append(count)
        Source.source_count +=1           
        Source.source_list.append(self)
        
    def add_pixel(self,location,count):
        self.pixel_coords.append(location)
        self.pixel_counts.append(count)
        

def mask_frame(mask,frame_cutoffs):

    for j in range(0,frame_y_top):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(( (dimensions[0]) -frame_y_bottom),dimensions[0]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
    
    mask = np.transpose(mask)
    
    for j in range(0,frame_x_left):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(( (dimensions[1]) -frame_x_right),dimensions[1]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    mask=np.transpose(mask)
    
    return mask

def source_detection(data,mask,background,dimensions,source_map,frame_cutoffs): 
    # Cyles through each pixel checking for sources
    count = 0
    for y in range(0,len(data)):
        for x in range(0,len(data[0])):            
            
            if mask[y][x] == True:            
                if data[y][x] > background:
                    count, source_map = source_handler(data,x,y,mask,background,dimensions,count,source_map)
    return count,mask,source_map
                
def source_handler(data,x,y,mask,background,dimensions,count,source_map):
    # Once a source is found, need to determine if pixel is part of a mulit-pixel distributed source
    x_low = x-1
    x_high = x+1
    y_high = y+1
    y_low = y-1

    if x_low < 0:
        x_low = 0
    if x_high >= dimensions[1]:
        x_high = dimensions[1] -1
    if y_low < 0:
        y_low = 0
    if y_high >= dimensions[0]:
        y_high = dimensions[0] -1                    

    source_added = False
    for j in range(y_low,y_high+1):         
        for i in range(x_low,x_high+1):
            
            if mask[j][i] == True:
            
                if source_map[j][i] > 0:
                    source_map[y][x] = source_map[j][i]
                    
                    if source_added == False:
                        index = int(source_map[j][i] -1)
                        Source.source_list[index].add_pixel((x,y), data[y][x])
                        source_added = True
                    
    if  source_map[y][x] == 0:
        count +=1
        s = Source(count,(x,y), data[y][x])

        source_map[y][x] = count
    return count, source_map

mask = mask_frame(mask, frame_cutoffs)
count, mask, source_map = source_detection(data,mask,background,dimensions,source_map,frame_cutoffs)





print("MASK: \n",mask)
print("SOURCE MAP: \n", source_map)
print("OBJECT COUNT:",count)