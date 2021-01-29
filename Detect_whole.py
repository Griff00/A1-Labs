# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 12:29:10 2021

@author: Jack
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pickle as pk

hdulist = fits.open("mosaic.fits")
real_data = hdulist[0].data
header = hdulist[0].header
inst_zero_point = header['MAGZPT']
inst_zero_point_err = header['MAGZRR']

generated_test_data = [[1,1,1,2,1,2,1,2,1,2,1,1],
        [1,1,1,2,1,30,30,2,1,2,40,1],
        [1,1,1,2,1,30,30,2,1,2,1,1],
        [1,2,1,2,1,2,1,2,1,2,1,1]]

#Cutting out the specified test region
#test_data = list(map(lambda x: x[660:1260], real_data[3960:4340]))
#data = test_data
data = real_data

#-------------------------------------------------------#
#----------------------HISTOGRAM------------------------#
#-------------------------------------------------------#
plot_histogram = False
if plot_histogram:
    #Change this to make the histogram more/less detailed
    histogram_no_of_bins = 500
    
    ''' 
    mins = []
    for row in data:
        mins.append(min(row))
    data_min = min(mins)
    
    maxes = []
    for row in data:
        maxes.append(max(row))
    data_max = max(maxes)
    '''
    
    #Use these pre-set values as the whole dataset was skewing the Gaussian
    data_min = 2600
    data_max = 3600
    
    #Set up the bins
    histogram_bin_size = (data_max - data_min)/histogram_no_of_bins
    histogram_bin_lowest_vals = list(range(data_min,data_max,int(histogram_bin_size)))[:histogram_no_of_bins]
    
    #Count the number of pixels in each bin
    histogram_counts = np.zeros(histogram_no_of_bins)
    for row in data:
        for i in row:
            bin_number = int(np.floor((i - data_min) / histogram_bin_size))
            if bin_number == histogram_no_of_bins:
                bin_number -= 1
            if bin_number < histogram_no_of_bins and not bin_number < 0 : 
                histogram_counts[bin_number] += 1
                
    plt.plot(histogram_bin_lowest_vals,histogram_counts, 'kx')
    plt.xlabel('Number of counts')
    plt.ylabel('Number of incidences')
    
    #Fit Gaussian to noise
    
    def Gaussian(x, mean, std, scaling):
        return scaling*np.e**(-(x-mean)**2/(2*std**2))
    
    #The actual fitting to the data - settles but errors are v high- most of the vals being considered are 0. 
    gaussian_params, pcov= spo.curve_fit(Gaussian,histogram_bin_lowest_vals,histogram_counts, p0=[3240,200,1e7])
    gaussian_params_errs = np.sqrt(np.diag(pcov))
    
    #Generating a plottable version 
    g_xvals = np.linspace(data_min, data_max, 1000)
    g_yvals = Gaussian(g_xvals, *gaussian_params)
    
    plt.plot(g_xvals,g_yvals, 'r-')
    plt.show()
    
    print("The parameters used to calculate this fit were a mean of %f \pm %f, a standard deviation of %f \pm %f, and a scaling factor of %f \pm %f"%(gaussian_params[0], gaussian_params_errs[0], gaussian_params[1], gaussian_params_errs[1] , gaussian_params[2],gaussian_params_errs[2]))


#-----------------------------------------------------------#
#------------------SOURCE DETECTION-------------------------#
#-----------------------------------------------------------#


class Source:
    #List of all sources (regardless if they are galaxies or not)
    source_list = []
    #j: Dictionary of [id:source_object]
    source_dict = {}
    #List of all objects identified as galaxies
    galaxy_list = []
    # (temporary) lists for debugging purposes
    segment_list = []
    removed_list = []
    allobs_list = []
    
    def __init__(self,identifier,location,count):                 
        self.num_pixels = 1
        self.identifier = identifier
        self.pixel_coords = [location]
        self.pixel_counts = [count] 
        self.pixel_count_errors = [np.sqrt(count)]
        self.background_contribution = [0,0]
        self.mag = np.inf
        self.mag_err = 0
        
        Source.source_list.append(self)
        Source.galaxy_list.append(self)
        Source.allobs_list.append(self)
        Source.source_dict[identifier] = self
        
    def add_pixel(self,location,count):
        self.num_pixels +=1
        self.pixel_coords.append(location)
        self.pixel_counts.append(count)
        self.pixel_count_errors.append(np.sqrt(count))
        
    def set_background_contribution(self,background_contribution, background_contribution_err):
        self.background_contribution = background_contribution
        self.background_contribution_err = background_contribution_err

    def remove_source(self,mask,galaxy_map):
        Source.galaxy_list.remove(self)
        Source.removed_list.append(self)
                            
        coords = self.get_pixel_coords()
        for i in range(0,len(coords)):
            galaxy_map[coords[i][1]][coords[i][0]] = 0
            mask[coords[i][1]][coords[i][0]] = False 
            
        return mask, galaxy_map

    def get_pixel_coords(self):
        return self.pixel_coords

    def get_num_pixels(self):
        return self.num_pixels
    
    def get_pixel_counts(self):
        return self.pixel_counts
   
    def find_magnitude(self, inst_zero_point, inst_zero_point_err):
        tot_counts = sum(self.pixel_counts) - self.background_contribution
        tot_counts_err = np.sqrt(sum(map(lambda x: x**2, self.pixel_count_errors)) + self.background_contribution_err**2)
        
        if tot_counts <= 0:
            print("tot",tot_counts,"num_pix",self.num_pixels)
            print("counts",self.pixel_counts,"back",self.background_contribution)
            tot_counts = 1
            tot_counts_err = 0
        
        self.mag = inst_zero_point -2.5*np.log(tot_counts)/np.log(10)
        self.mag_err = np.sqrt(inst_zero_point_err**2 + (2.5/np.log(10)/tot_counts*tot_counts_err)**2)
        
        return self.mag, self.mag_err

def mask_frame(mask,frame_cutoffs):
    for j in range(0,frame_cutoffs["yt"]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(((len(data))-frame_cutoffs["yb"]),len(data)):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
    
    mask = np.transpose(mask)
    
    for j in range(0,frame_cutoffs["xl"]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(( (len(data[0])) -frame_cutoffs["xr"]),len(data[0])):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    mask=np.transpose(mask)
    
    return mask

def mask_box(mask, box_x1, box_x2, box_y1, box_y2):
    for i in range(box_x1, box_x2):
        for j in range(box_y1,box_y2):
            mask[j][i] = False
    return mask

def calculate_background_contribution(source, data,mask, global_background,source_map, aperture_radius=12):
    #Calculates the total background contribution to a given source
    n_pixels_considered = 0
    counts_considered = []
    coords = source.get_pixel_coords()
    
    #calculate average centres for the surce (where aperture will be centred)
    x_avg_source = sum(map(lambda x: x[0], coords))/len(coords)
    y_avg_source = sum(map(lambda x: x[1], coords))/len(coords)
    
    x_low = int(x_avg_source - (aperture_radius + 5))
    x_high = int(x_avg_source + (aperture_radius + 5))
    y_high = int(y_avg_source + (aperture_radius + 5))
    y_low = int(y_avg_source - (aperture_radius + 5))
    
    if x_low < 0:
        x_low = 0
    if x_high >= len(data[0]):
        x_high = len(data[0]) -1
    if y_low < 0:
        y_low = 0
    if y_high >= len(data):
        y_high = len(data) -1
    
    #find + sum over the aperture
    for y in range(y_low,y_high):
        for x in range(x_low,x_high):            
            if (y-y_avg_source)**2+(x-x_avg_source)**2 < aperture_radius**2:
                if not (x,y) in coords and mask[y][x] == True and source_map[y][x] == 0:
                    counts_considered.append(data[y][x])
                    n_pixels_considered += 1
                    
    if n_pixels_considered != 0:
        local_bg_emmission = sum(counts_considered)/n_pixels_considered
        local_emmission_err = np.std(counts_considered)
        total_contribution_to_source = local_bg_emmission*len(coords)
    else:
        total_contribution_to_source = global_background*len(coords)
        local_emmission_err = np.sqrt(global_background)
        
    
    return total_contribution_to_source, local_emmission_err
 

   

def source_detection(data,mask,background,source_map,frame_cutoffs): 
    # Cyles through each pixel checking for sources
    source_id = 0
    for y in range(0,len(data)):
        for x in range(0,len(data[0])):
            if mask[y][x] == True:            
                if data[y][x] > background:
                    #If source pixel found --> call source_handler
                    source_id, source_map = source_handler(data,x,y,mask,source_id,source_map)
    return source_id,mask,source_map

def source_handler(data,x,y,mask,source_id,source_map):
    # Once a source is found, need to determine if pixel is part of a mulit-pixel distributed source
    #Check neighbour pixels for already catalogued sources
    x_low = x-1
    x_high = x+1
    y_high = y+1
    y_low = y-1

    #Dealing with edge cases to avoid 'out of range' errors
    if x_low < 0:
        x_low = 0
    if x_high >= len(data[0]):
        x_high = len(data[0]) -1
    if y_low < 0:
        y_low = 0
    if y_high >= len(data):
        y_high = len(data) -1              

    for j in range(y_low,y_high+1):         
        for i in range(x_low,x_high+1):
            #Check if pixel is masked
            if mask[j][i] == True:
                #If there is neighbour source
                if source_map[j][i] > 0: 
                    #Add current pixel to neighbour source
                    if source_map[y][x] == 0:
                        source_map[y][x] = source_map[j][i]
                        Source.source_dict[source_map[y][x]].add_pixel((x,y), data[y][x])
                    #Pixel has at least 2 neighbour sources --> need to combine these neighbour segments into 1 source 
                    elif source_map[y][x] != source_map[j][i]:
                        #Identify 'true' source, and the 'segment' to be added (true source has earliest id)
                        if source_map[y][x] > source_map[j][i]:
                            true_id = source_map[j][i]
                            segment_id = source_map[y][x]
                        else:
                            true_id = source_map[y][x]
                            segment_id = source_map[j][i]
                    
                        #Retrieve source objects & remove segment source from source_list
                        true_source = Source.source_dict[true_id]
                        segment_source = Source.source_dict[segment_id]
                        Source.source_list.remove(segment_source)
                        Source.segment_list.append(segment_source)
                        
                        coords = segment_source.get_pixel_coords()
                        pixel_counts = segment_source.get_pixel_counts()
                        
                        #Add all segment source pixel data to true source
                        for i in range(0,len(coords)):
                            true_source.add_pixel(coords[i],pixel_counts[i])
                            source_map[coords[i][1]][coords[i][0]] = true_id
                        
    # No neighbours --> new source              
    if  source_map[y][x] == 0:
        source_id +=1
        Source(source_id,(x,y), data[y][x])
        source_map[y][x] = source_id
    
    return source_id, source_map

def contam_removal(mask, source_map,lower_pix_limit, upper_pix_limit):
    # Remove contamination (e.g. single broken pixels, or significantly large artifacts)
    galaxy_map = source_map.copy()
    galaxy_list_orig = Source.galaxy_list.copy()
    
    for ob in galaxy_list_orig:        
        if ob.get_num_pixels() <= lower_pix_limit or ob.get_num_pixels() >= upper_pix_limit:
            #Remove source
            mask, galaxy_map = ob.remove_source(mask,galaxy_map)
                          
    return galaxy_map, mask
   

#-------------------BUILDING THE MASK---------------------------------#

mask = np.ones((len(data),len(data[0])), dtype=bool)
#For real data
#frame_cutoffs = {"xl":150,"xr":2450,"yt":4500,"yb":150} 

#Values are the thickness (in number of pixels) to cut off (not coordinate positions)
frame_cutoffs = {"xl":150,"xr":120,"yt":111,"yb":150} 
mask = mask_frame(mask, frame_cutoffs)

#areas with high noise in the corners
mask = mask_box(mask,0,410,0,200)
mask = mask_box(mask,0,200,0,410)
mask = mask_box(mask,2100,len(mask),4400,len(mask[0]))
mask = mask_box(mask,2350,len(mask),4100,len(mask[0]))

#bright central star (and bleeding from it)

# Original mask values:
# mask = mask_box(mask,1200,1650,2950,3450)

#Expanded mask values:
mask = mask_box(mask,1000,1800,2700,3700)
mask = mask_box(mask,1420,1460,0,len(mask[0]))
mask = mask_box(mask,1100,1650,0,500)

#--------------BACKGROUND THRESHOLD-----------------------------#

# Debug / testing threshold:



# Real threshold
sigma_away = 4
# From gaussian fit ---> so don't need to keep running gaussian
mean = 3418.2
std = 12.2

background = mean + sigma_away*std 
#background = gaussian_params[0] + sigma_away*abs(gaussian_params[1])

#j: lower & upper pixel number limits for removing contamination
lower_pix_limit = 3
upper_pix_limit = 1000

#---------------------------------------------------------------#

source_map = np.zeros((len(data),len(data[0])))
count, mask, source_map = source_detection(data,mask,background,source_map,frame_cutoffs)
galaxy_map, mask = contam_removal(mask, source_map,lower_pix_limit, upper_pix_limit)

iteration_list = Source.galaxy_list.copy()
for ob in iteration_list:
    background_contribution, bg_error = calculate_background_contribution(ob, data,mask,background,source_map, 12)
    ob.set_background_contribution(background_contribution, bg_error)
    
    mag = ob.find_magnitude(inst_zero_point, inst_zero_point_err)


#print("MASK: \n",mask)
#print("SOURCE MAP: \n", source_map)
#print("OBJECT COUNT:",count)
galax_count = len(Source.galaxy_list)
#print("GALAXY MAP: \n", galaxy_map)
print("GALAXY COUNT:",galax_count)
print("Sigma Away (used to set threshold)=",sigma_away)
print("THESHOLD =",background, "(counts per pixel)")
print("Minimum source pixel number",lower_pix_limit)
print("Maximum source pixel number",upper_pix_limit)
num_pix_unmasked =0
for i in range(0,len(mask)):
    for j in range(0,len(mask[0])):
        if mask[i][j] == True:
            num_pix_unmasked +=1
print("Number of (unmasked) pixels searched over =",num_pix_unmasked)


#---------Plot of galaxies (for debugging)-----------#
x = []
y = []
s = []
for ob in Source.galaxy_list:
    coords = ob.get_pixel_coords()
    size = ob.get_num_pixels()
    x_avg_source = sum(map(lambda x: x[0], coords))/len(coords)
    y_avg_source = sum(map(lambda x: x[1], coords))/len(coords)
    x.append(x_avg_source)
    y.append(y_avg_source)
    s.append(size)
    
plt.figure()
plt.scatter(x,y,s)
plt.xlabel('x pixel coordinates')
plt.ylabel('y pixel coordinates')
plt.title('Visual Representation (for Debugging) of identified Galaxies within image subsection. \n Galaxy coordinates and ~size shown.')
plt.show()



mags = []
mags_err = []
for ob in Source.galaxy_list:
    mags.append(ob.mag)
    mags_err.append(ob.mag_err)    
out_data = [mags,mags_err]
np.savetxt('galaxy_list.csv',out_data, delimiter=",")