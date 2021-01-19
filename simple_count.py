
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo


hdulist = fits.open("mosaic.fits")

#Test data
real_data = hdulist[0].data
generated_test_data = [[1,1,1,2,1,2,1,2,1,2,1,1],
        [1,1,1,2,1,30,30,2,1,2,40,1],
        [1,1,1,2,1,30,30,2,1,2,1,1],
        [1,2,1,2,1,2,1,2,1,2,1,1]]

header = hdulist[0].header

instrument_zero_point = header['MAGZPT']
instrument_zero_point_err = header['MAGZRR']

#Cutting out the specified test region
test_data = list(map(lambda x: x[660:1260], real_data[3960:4340]))

data = test_data

#----------------------HISTOGRAM------------------------#

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
histogram_no_of_bins = 100

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


background = gaussian_params[0] + gaussian_params[1] #Only an initial estimate

mask = np.ones((len(data), len(data[0])), dtype=bool)
source_map = np.zeros((len(data), len(data[0])))

#g: made below into a dictionary, so it keeps meaningful names while storing the same info
frame_cutoffs = {"xl":150,"xr":2450,"yt":4500,"yb":150} #for the full real data


#j: lower & upper pixel number limits for removing contamination
lower_pix_limit = 1
upper_pix_limit = 50

class Source:
    #g: not sure if this is the right place for this- we can make a list of the sources outside their class.
    #   would take plenty of restructuring so will not do it now, but smth to think abt
    source_list = []
    galaxy_list = []
    
    def __init__(self,identifier,location,count):
                 
        self.num_pixels = 1
        self.identifier = identifier
        self.pixel_coords = [location]
        self.pixel_counts = [count] 
        Source.source_list.append(self)
        Source.galaxy_list.append(self)
        
    def add_pixel(self,location,count):
        self.num_pixels +=1
        self.pixel_coords.append(location)
        self.pixel_counts.append(count)
        
    def get_pixel_coords(self):
        return self.pixel_coords

    def get_num_pixels(self):
        return self.num_pixels
    
    def get_magnitude(self, inst_zero_point):
        tot_counts = sum(self.pixel_counts)
        mag = inst_zero_point -2.5*np.log(tot_counts)/np.log(10)
        return mag

def mask_frame(mask,frame_cutoffs):
    #g: changed this so it was actually using the values passed in, instead of prestored values
    for j in range(0,frame_cutoffs["yt"]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(((len(mask))-frame_cutoffs["yb"]),len(mask)):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
    
    mask = np.transpose(mask)
    
    for j in range(0,frame_cutoffs["xl"]):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    for j in range(((len(mask)) -frame_cutoffs["xr"]),len(mask)):
        for i in range(0,len(mask[0])):
            mask[j][i] = False
            
    mask=np.transpose(mask)
    
    return mask

def mask_box(mask, box_x1, box_x2, box_y1, box_y2):
    for i in range(box_x1, box_x2):
        for j in range(box_y1,box_y2):
            mask[i][j] = False
    return mask

def calculate_background_contribution(source, data, mask, aperture_radius=12): #update so checks against mask
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
                if not (x,y) in coords and mask[x][y] == True:
                    total_counts += data[y][x]
                    n_pixels_considered += 1
                    
    local_bg_emmission = total_counts/n_pixels_considered
    total_contribution_to_source = local_bg_emmission*len(coords)
    
    return total_contribution_to_source
    

def source_detection(data,mask,background,source_map,frame_cutoffs): 
    # Cyles through each pixel checking for sources
    count = 0
    for y in range(0,len(data)):
        for x in range(0,len(data[0])):
            if mask[y][x] == True:      
                if data[y][x] > background:
                    count, source_map = source_handler(data,x,y,mask,count,source_map)
    return count,mask,source_map
                
def source_handler(data,x,y,mask,count,source_map):
    # Once a source is found, need to determine if pixel is part of a mulit-pixel distributed source
    x_low = x-1
    x_high = x+1
    y_high = y+1
    y_low = y-1

    if x_low < 0:
        x_low = 0
    if x_high >= len(data[0]):
        x_high = len(data[0]) -1
    if y_low < 0:
        y_low = 0
    if y_high >= len(data):
        y_high = len(data) -1              

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
        Source(count,(x,y), data[y][x])

        source_map[y][x] = count
    return count, source_map

def contam_removal(mask, source_map,lower_pix_limit, upper_pix_limit):
    #print("running")
    galaxy_map = source_map.copy()
    for ob in Source.galaxy_list:
        if ob.get_num_pixels() <= lower_pix_limit or ob.get_num_pixels() >= upper_pix_limit:
            Source.galaxy_list.remove(ob)
            
            for j in range (0,len(galaxy_map)):
                for i in range(0,len(galaxy_map[0])):
                    if galaxy_map[j][i] == ob.identifier:
                        galaxy_map[j][i] = 0
                        mask[j][i] = False
                        #print("here")
                        
    return galaxy_map, mask

mask = mask_frame(mask, frame_cutoffs)

#areas with high noise in the corners
mask = mask_box(mask,0,410,0,200)
mask = mask_box(mask,0,200,0,410)
mask = mask_box(mask,2100,len(mask),4400,len(mask[0]))
mask = mask_box(mask,2350,len(mask),4100,len(mask[0]))

#bright central star (and bleeding from it)
mask = mask_box(mask,1200,1650,2950,3450)
mask = mask_box(mask,1420,1460,0,len(mask[0]))
mask = mask_box(mask,1100,1650,0,500)


count, mask, source_map = source_detection(data,mask,background,source_map,frame_cutoffs)
galaxy_map, mask = contam_removal(mask, source_map,lower_pix_limit, upper_pix_limit)


print("MASK: \n",mask)
print("SOURCE MAP: \n", source_map)
print("OBJECT COUNT:",count)
print("GALAXY MAP: \n", galaxy_map)