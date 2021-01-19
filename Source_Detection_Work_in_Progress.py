
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo


hdulist = fits.open("mosaic.fits")
header = hdulist[0].header
instrument_zero_point = header['MAGZPT']
instrument_zero_point_err = header['MAGZRR']
data_original = hdulist[0].data


# data = [[1,1,1,2,1,2,1,2,1,2,1,1],
#         [1,1,1,2,1,30,30,2,1,2,40,1],
#         [1,1,1,2,1,30,30,2,1,2,1,1],
#         [1,2,1,2,1,2,1,2,1,2,1,1]]


#--------------SMALL PORTION OF IMAGE----------------#
dimensions_orig = (len(data_original),len(data_original[0]))
data1 = data_original[dimensions_orig[0] - 4500 : dimensions_orig[0] - 3000]
data1 = np.transpose(data1)
data2 = data1[500:2000]
data2 =np.transpose(data2)
#j: data is smaller portion of data original
data = data2.copy()

dimensions = (len(data),len(data[0]))


#-------------------------------------------------------#
#----------------------HISTOGRAM------------------------#
#-------------------------------------------------------#


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




#-----------------------------------------------------------#
#------------------SOURCE DETECTION-------------------------#
#-----------------------------------------------------------#


mean = np.mean(data)
std = np.std(data)

background = gaussian_params[0] + 3*gaussian_params[1]
#background = mean + std #Just tempory (for testing)

mask = np.ones(dimensions, dtype=bool)
source_map = np.zeros(dimensions)

#------------------MASK PARAMETERS-------------------------#

frame_x_left = 1
frame_x_right = 1
frame_y_top = 1
frame_y_bottom = 1
#g: made below into a dictionary, so it keeps meaningful names while storing the same info
frame_cutoffs = {"xl":frame_x_left,"xr":frame_x_right,"yt":frame_y_top,"yb":frame_y_bottom}

#j: lower & upper pixel number limits for removing contamination
lower_pix_limit = 1
upper_pix_limit = 100000000

#--------------------------------------------------------#

class Source:
    #g: not sure if this is the right place for this- we can make a list of the sources outside their class.
    #   would take plenty of restructuring so will not do it now, but smth to think abt
    source_list = []
    #j: Dictionary of id:source_object
    source_dict = {}
    galaxy_list = []
    # lists for debugging purposes
    segment_list = []
    removed_list = []
    allobs_list = []
    
    def __init__(self,identifier,location,count):
                 
        self.num_pixels = 1
        self.identifier = identifier
        self.pixel_coords = [location]
        self.pixel_counts = [count] 
        Source.source_list.append(self)
        Source.galaxy_list.append(self)
        Source.allobs_list.append(self)
        Source.source_dict[identifier] = self
        
    def add_pixel(self,location,count):
        self.num_pixels +=1
        self.pixel_coords.append(location)
        self.pixel_counts.append(count)
        
    def get_pixel_coords(self):
        return self.pixel_coords

    def get_num_pixels(self):
        return self.num_pixels
    
    def get_pixel_counts(self):
        return self.pixel_counts

    def get_id(self):
        return self.identifier
        



def mask_frame(mask,frame_cutoffs):
    #g: changed this so it was actually using the values passed in, instead of prestored values
    for j in range(0,frame_cutoffs["yt"]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(((dimensions[0])-frame_cutoffs["yb"]),dimensions[0]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
    
    mask = np.transpose(mask)
    
    for j in range(0,frame_cutoffs["xl"]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(( (dimensions[1]) -frame_cutoffs["xr"]),dimensions[1]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    mask=np.transpose(mask)
    
    return mask



def calculate_background_contribution(source, data, aperture_radius=12):
    #Calculates the total background contribution to a given source
    n_pixels_considered = 0
    total_counts = 0
    coords = source.get_pixel_coords()
    
    #calculate average centres for the surce (where aperture will be centred)
    x_avg_source = sum(map(lambda x: x[0], coords))/len(coords)
    y_avg_source = sum(map(lambda x: x[1], coords))/len(coords)
    
    #find + sum over the aperture
    for y in range(0,len(data)):
        for x in range(0,len(data[0])):
            if (y-y_avg_source)**2+(x-x_avg_source)**2 < aperture_radius**2:
                if not (x,y) in coords:
                    total_counts += data[y][x]
                    n_pixels_considered += 1
                    
    local_bg_emmission = total_counts/n_pixels_considered
    total_contribution_to_source = local_bg_emmission*len(coords)
    
    return total_contribution_to_source
    



def source_detection(data,mask,background,dimensions,source_map,frame_cutoffs): 
    # Cyles through each pixel checking for sources
    source_id = 0
    for y in range(0,len(data)):
        for x in range(0,len(data[0])):
            if mask[y][x] == True:            
                if data[y][x] > background:
                    #If source pixel found --> call source_handler
                    source_id, source_map = source_handler(data,x,y,mask,dimensions,source_id,source_map)
    return source_id,mask,source_map


               
def source_handler(data,x,y,mask,dimensions,source_id,source_map):
    # Once a source is found, need to determine if pixel is part of a mulit-pixel distributed source
    #Check neighbour pixels for already catalogued sources
    x_low = x-1
    x_high = x+1
    y_high = y+1
    y_low = y-1


    #Dealing with edge cases to avoid 'out of range' errors
    if x_low < 0:
        x_low = 0
    if x_high >= dimensions[1]:
        x_high = dimensions[1] -1
    if y_low < 0:
        y_low = 0
    if y_high >= dimensions[0]:
        y_high = dimensions[0] -1              

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
                        
                        if segment_id == 691:
                            print("here")
                        
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
        
        if source_id == 691:
            print("create")
        
        Source(source_id,(x,y), data[y][x])

        source_map[y][x] = source_id
    return source_id, source_map

def contam_removal(mask, source_map,lower_pix_limit, upper_pix_limit):
    # Remove contamination (e.g. single broken pixels, or significantly large artifacts)
    galaxy_map = source_map.copy()
    
    print("bool",Source.source_dict[691] in Source.galaxy_list)
    
    for ob in Source.galaxy_list:
        
        if ob.get_id() > 400 and ob.get_id() < 700:
            print("ob id",ob.get_id())
        
        if ob.get_num_pixels() <= lower_pix_limit or ob.get_num_pixels() >= upper_pix_limit:
            Source.galaxy_list.remove(ob)
            Source.removed_list.append(ob)
                   
            if ob.get_id() == 691:
                print("found")
             
            coords = ob.get_pixel_coords()
            for i in range(0,len(coords)):
                #print("coords",coords)
                
                galaxy_map[coords[i][1]][coords[i][0]] = 0
                mask[coords[i][1]][coords[i][0]] = False
    print("bool2",Source.source_dict[691] in Source.galaxy_list)                             
    return galaxy_map, mask

mask = mask_frame(mask, frame_cutoffs)
count, mask, source_map = source_detection(data,mask,background,dimensions,source_map,frame_cutoffs)
galaxy_map, mask = contam_removal(mask, source_map,lower_pix_limit, upper_pix_limit)


print("MASK: \n",mask)
print("SOURCE MAP: \n", source_map)
print("OBJECT COUNT:",count)
galax_count = len(Source.galaxy_list)
print("GALAXY MAP: \n", galaxy_map)
print("GALAXY COUNT:",galax_count)