"""
Created on Sun Dec 27 18:31:22 2020
"""
#%%
#Validation and graphs created for the astro imaging code

import matplotlib.pyplot as plt

hdulist = fits.open("/Users/philippsonnenschein/Desktop/Laboratory/Astro_Imaging/A1_mosaic.fits")
#hdulist = fits.open('A1_mosaic.fits')

MAGZPT = hdulist[0].header[157] 
MAGZRR = hdulist[0].header[158] 
print(MAGZPT,"+/-",MAGZRR)

linear_data = raw_data.flatten()    #Length 11850270
raw_data = hdulist[0].data

y_len = 4611
x_len = 2570

hdulist.close("/Users/philippsonnenschein/Desktop/Laboratory/Astro_Imaging/A1_mosaic.fits")
#hdulist.close('A1_mosaic.fits') 

#%%
#Data Validation

#Plotting the distribution of values
#Careful: This takes a lot of time and is very slow

plt.hist(filt_data, bins=100)
plt.xlabel('Pixel value')
plt.ylabel('intensity')
plt.show()

#%%
#Background Error distribution

filt_linear_data = filt_data.flatten()
error_data=[]
for i in filt_linear_data:      #selecting the relevant data
    if 3300<i<3700:
        error_data.append(i)

standard_dev=np.std(error_data)
mean=np.mean(error_data)
#standard_dev = 12.957414307332735
#mean = 3419.8636867360174
print("Standard deviation =",standard_dev)
print("Mean =",mean)

def normal_distribution(mu,std,x):
    r = 2 * np.pi
    n = 1/(std*np.sqrt(r))
    exponential = np.exp((-(x-mu)**2)/(2*std**2))
    normal_function = n * exponential
    return normal_function

#%%
#Error Distribution according to determined std and mean

print("Standard Deviation =",standard_dev)
print("Mean =",mean)

Gaussian = normal_distribution(mean,standard_dev,error_data)
Gaussian_scaled = Gaussian * 14335000 * (288700/441600)               #Scaling is arbitrary to fit hist magnitude
plt.plot(error_data, Gaussian_scaled,'rx', markersize=2)
plt.hist(error_data, bins=500)
plt.xlabel('Pixel Value')
plt.ylabel('Number of Occurences')
plt.show()   

#%%
#Background noise confidence value evaluation

set_threshhold = 3500

from scipy import integrate

x1 = np.linspace(0,6000,100000)
Gaussian1 = normal_distribution(mean,standard_dev,x1)
x2 = np.linspace(0,set_threshhold,100000)
Gaussian2 = normal_distribution(mean,standard_dev,x2)


area_full_gaussian = integrate.simps(Gaussian1,x1)
area_covered = integrate.simps(Gaussian2,x2)
print("There is a ",(1 - (area_covered/area_full_gaussian))*100,"Percent chance that the detected galaxy is actually background noise")


print((1-(area_covered/area_full_gaussian))*variable_num_galaxies,"Galaxies can be attributed to background noise")

#5-sigma confidence value implies an error of 0.00006%

#%%
#Fixed Aperture Validation

#Plot of the locations of the galaxies

size = len(fixed_catalogue)
fixed_x_coordinates = np.zeros(size)
fixed_y_coordinates = np.zeros(size)
for i in range(0,size):
    fixed_x_coordinates[i] = fixed_catalogue[i][2]
    fixed_y_coordinates[i] = fixed_catalogue[i][3]
    
plt.figure()
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Galaxies, Fixed Aperture Method')
plt.scatter(fixed_x_coordinates,fixed_y_coordinates,marker='x')
plt.show()

#%%
#Plot of the locations of the galaxies with sizes weighted by the determined radius

size = len(fixed_catalogue)
fixed_x_coordinates = np.zeros(size)
fixed_y_coordinates = np.zeros(size)
fixed_radius = np.zeros(size)
for i in range(0,size):
    fixed_x_coordinates[i] = fixed_catalogue[i][2]
    fixed_y_coordinates[i] = fixed_catalogue[i][3]
    fixed_radius[i] = fixed_catalogue[i][6]
#plt.figure()
plt.axis('equal')
plt.xlim(0,2600)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Galaxies, Fixed Aperture Method')
for i in range(size):
    plt.scatter(fixed_x_coordinates[i],fixed_y_coordinates[i],c='b',marker='o',s=int(fixed_radius[i]**2 /8))
plt.show()

#%%
#Plot the distribution of the local backgrounds
#compare to the gaussian distribution from before

size = len(fixed_catalogue)

max_background = 0
for i in range(size-1):
    if fixed_catalogue[i][7] > max_background:
        max_background = fixed_catalogue[i][7]
print(max_background)

min_background = 4000
for i in range(size-1):
    if fixed_catalogue[i][7] < min_background:
        if fixed_catalogue[i][7] != 0:
            min_background = fixed_catalogue[i][7]
print(min_background)

#Max is 3481
#Min is 3408
#Difference is 73


Background = np.linspace(int(np.floor(min_background)),int(np.ceil(max_background)),int((np.ceil(max_background))-np.floor(min_background)+1),endpoint='True')

Mag_backgound = np.zeros(int((np.ceil(max_background))-np.floor(min_background)+1))
for i in range(1,size-1):
    back_g = int(fixed_catalogue[i][7])-int(np.floor(min_background))
    Mag_backgound[back_g] = Mag_backgound[back_g] + 1

plt.xlim(3350,3550)
plt.scatter(Background,Mag_backgound,marker='x')


#measured_mean = np.std(y_values4)
#measured_standard_dev = np.mean(y_values4)
mean = 3419
standard_dev = 6
x_values = np.linspace(3000,3700,100000,dtype=int)

Gaussian = normal_distribution(mean, standard_dev, x_values)
Gaussian_scaled = Gaussian*1900
plt.plot(x_values, Gaussian_scaled,'gx', markersize=2)
plt.xlabel('Pixel value')
plt.ylabel('Number of Occurences')
plt.title('Distribution of the Local Background')
plt.show()     

#%%
#Variable Aperture Validation

#Plot of the locations of the galaxies

size = len(variable_catalogue)
variable_x_coordinates = np.zeros(size)
variable_y_coordinates = np.zeros(size)
for i in range(0,size):
    variable_x_coordinates[i] = variable_catalogue[i][2]
    variable_y_coordinates[i] = variable_catalogue[i][3]

plt.figure()
plt.axis('equal')
plt.xlim(0,2600)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Galaxies, Variable Aperture Method')
plt.scatter(variable_x_coordinates,variable_y_coordinates,marker='x')
plt.show()

#%%
#Plot of the locations of the galaxies with sizes weighted by the determined radius

size = len(variable_catalogue)
variable_x_coordinates = np.zeros(size)
variable_y_coordinates = np.zeros(size)
variable_radius = np.zeros(size)
for i in range(0,size):
    variable_x_coordinates[i] = variable_catalogue[i][2]
    variable_y_coordinates[i] = variable_catalogue[i][3]
    variable_radius[i] = variable_catalogue[i][6]
#plt.figure()
plt.axis('equal')
plt.xlim(0,2600)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Galaxies, Variable Aperture Method')
for i in range(size):
    plt.scatter(variable_x_coordinates[i],variable_y_coordinates[i],c='b',marker='o',s=int(variable_radius[i]**2 /8))
plt.show()

#%%
#Plot of the distribution of radii

size = len(variable_catalogue)
Rad = np.linspace(0,106,107) #the largest radius is 106 pixels
Mag_rad = np.zeros(107)
for i in range(0,size-1):
    radius = int(variable_catalogue[i][6])
    Mag_rad[radius] = Mag_rad[radius] + 1
    
plt.xlim(0,20)
plt.xlabel('Radius (Pixels)')
plt.ylabel('Number of Occurences')
plt.title('Distribution of the Radii')
plt.xticks([0,4,8,12,16,20])
#plt.yticks([200,250,300,350])
plt.scatter(Rad,Mag_rad,marker='x')
plt.show()

#%%
#Number of galaxies with a radius larger than 12 pixels

num_large_rad_galaxies = int(np.sum(Mag_rad[13:]))
print("There are ",num_large_rad_galaxies,"Galaxies with more than 12 pixels")
print("These pixels strongly increase the counted number of galaxies counted using the fixed aperture method")

#%%
#Plot the distribution of the local backgrounds
#compare to the gaussian distribution from before

size = len(variable_catalogue)

max_background = 0
for i in range(size-1):
    if variable_catalogue[i][7] > max_background:
        max_background = variable_catalogue[i][7]
print(max_background)

min_background = 4000
for i in range(size-1):
    if variable_catalogue[i][7] < min_background:
        if variable_catalogue[i][7] != 0:
            min_background = variable_catalogue[i][7]
print(min_background)

#Max is 3481
#Min is 3408
#Difference is 73


Background = np.linspace(int(np.floor(min_background)),int(np.ceil(max_background)),int((np.ceil(max_background))-np.floor(min_background)+1),endpoint='True')

Mag_backgound = np.zeros(int((np.ceil(max_background))-np.floor(min_background)+1))
for i in range(1,size-1):
    back_g = int(variable_catalogue[i][7])-int(np.floor(min_background))
    Mag_backgound[back_g] = Mag_backgound[back_g] + 1

plt.xlim(3350,3550)
plt.scatter(Background,Mag_backgound,marker='x')


#measured_mean = np.std(y_values4)
#measured_standard_dev = np.mean(y_values4)
mean = 3436
standard_dev = 10
x_values = np.linspace(3000,3700,100000,dtype=int)

Gaussian = normal_distribution(mean, standard_dev, x_values)
Gaussian_scaled = Gaussian*2700
plt.plot(x_values, Gaussian_scaled,'gx', markersize=2)
plt.xlabel('Pixel value')
plt.ylabel('Number of Occurences')
plt.title('Distribution of the local background')
plt.show()     

#%% Check for eliptical galaxies
#Assume a threshold of 4 pixels over which the galaxy can be classified as elliptical

size = len(variable_catalogue)
initial_count = 0       #Number of galaxies where the x diameter != y diameter
final_count = 0         #Number of galaxies where the difference is > 4 pixels
diam_dif = np.array([]) #Array of differences
for i in range(size):
    if variable_catalogue[i][8] != variable_catalogue[i][9]:
        difference = int(np.absolute(variable_catalogue[i][8] - variable_catalogue[i][9]))
        diam_dif = np.append(diam_dif,[difference])
        initial_count += 1
        if difference >= 4:
            final_count += 1
            
#print(initial_count)
print("Number of ellitical galaxies = ",final_count)

#Plot of the Different in Measured Diameters
max_dif = int(np.max(diam_dif))

diam_difference = np.linspace(1,max_dif,max_dif+1)
Mag_diam_difference = np.zeros(max_dif+1)
for i in range(len(diam_dif)):
    difference = int(diam_dif[i])
    Mag_diam_difference[difference] = Mag_diam_difference[difference] + 1
    
    
plt.xlim(1,max_dif)
plt.xlabel('Difference in Measured Diameter (pixels)')
plt.ylabel('Number of Occurrences')
plt.title('Plot of the Difference in Measured Diameters')
plt.scatter(diam_difference,Mag_diam_difference,marker='x')
plt.show()

#%%
#Number Count plot Variable Aperture

def linear_func(x,slope,yintercept):
    return (slope*x) + yintercept

MAGZPT = hdulist[0].header[157] 
MAGZRR = hdulist[0].header[158] 

variable_magnitude = MAGZPT-(2.5*(np.log10(variable_catalogue[:,4])))        #chnage from pixel number to magnitude
variable_unique_mag, variable_counts = np.unique(variable_magnitude, return_counts=True)
variable_cumulative = np.cumsum(variable_counts)
variable_log_number = np.log10(variable_cumulative)

#Error analysis
error_variable_magnitude = (2.5*(np.log10(0.1*variable_catalogue[:,4])))
error_variable_unique_mag, error_variable_counts = np.unique(error_variable_magnitude, return_counts=True)
error_variable_cumulative = np.cumsum(error_variable_counts)
error_variable_log_number = np.log10(error_variable_cumulative)

#The first 3 entries are wrong, could delete here or at the beginning when filtering data
#initial_index = [0,1,2]
initial_index = []
final_variable_unique_mag = np.delete(variable_unique_mag,initial_index)
final_variable_log_number = np.delete(variable_log_number,initial_index)

#Modified without the final entries where the count flattens out for an improved linear fit
final_index = np.arange(1700,2219)
mod_variable_unique_mag = np.delete(final_variable_unique_mag,final_index)
mod_variable_log_number = np.delete(final_variable_log_number,final_index)

theoretical_slope = 0.6
theoretical_yintercept = -5

print('Theoretical linear function:')
print('Slope = %.3e' %(theoretical_slope))# Note the %.3e forces the values to be printed in scientific notation with 3 decimal places.
print('Intercept = %.3e' %(theoretical_yintercept))
print()

yerror1 = 1.0/error_variable_log_number
errors1 = 1.0/mod_variable_log_number

fit1,cov1 = np.polyfit(mod_variable_unique_mag,mod_variable_log_number,1,w=1/errors1,cov=True)

measured_slope1 = fit1[0]
measured_yintercept1 = fit1[1]
sig_0_1 = np.sqrt(cov1[0,0]) #The uncertainty in the measured slope
sig_1_1 = np.sqrt(cov1[1,1]) #The uncertainty in the measured intercept

print('Fited linear function:')
print('Slope = %.3e +/- %.3e' %(measured_slope1,sig_0_1))# Note the %.3e forces the values to be printed in scientific notation with 3 decimal places.
print('Intercept = %.3e +/- %.3e' %(measured_yintercept1,sig_1_1))

plt.figure()
plt.errorbar(final_variable_unique_mag,final_variable_log_number,yerr=yerror1,fmt='x',capsize=2,label='Data')
plt.plot(mod_variable_unique_mag, linear_func(mod_variable_unique_mag,measured_slope1,measured_yintercept1),label='Linear Fit')
plt.plot(mod_variable_unique_mag, linear_func(mod_variable_unique_mag,theoretical_slope,theoretical_yintercept),label='Theoretical Model')

plt.xlim(7,21)
plt.ylim(-2,7)
plt.xlabel('Magnitude')
plt.ylabel('Log10(Count)')
plt.title('Number Count Plot using the Variable Aperture Method')
#plt.xticks([200,250,300,350])
#plt.yticks([200,250,300,350])
plt.legend(loc='best')
plt.savefig('Variable_Numbercount_Plot')
plt.show()



#%%
#Number Count plot Fixed Aperture

def linear_func(x,slope,yintercept):
    return (slope*x) + yintercept

MAGZPT = hdulist[0].header[157] 
MAGZRR = hdulist[0].header[158] 

fixed_magnitude = MAGZPT-(2.5*(np.log10(fixed_catalogue[:,4])))        #chnage from pixel number to magnitude
fixed_unique_mag, fixed_counts = np.unique(fixed_magnitude, return_counts=True)
fixed_cumulative = np.cumsum(fixed_counts)
fixed_log_number = np.log10(fixed_cumulative)

#Error analysis
error_fixed_magnitude = (2.5*(np.log10(0.1*fixed_catalogue[:,4])))
error_fixed_unique_mag, error_fixed_counts = np.unique(error_fixed_magnitude, return_counts=True)
error_fixed_cumulative = np.cumsum(error_fixed_counts)
error_fixed_log_number = np.log10(error_fixed_cumulative)


#The first 3 entries are wrong, could delete here or at the beginning when filtering data
#initial_index = [0,1,2]
initial_index = []
final_fixed_unique_mag = np.delete(fixed_unique_mag,initial_index)
final_fixed_log_number = np.delete(fixed_log_number,initial_index)

#Modified without the final entries where the count flattens out for an improved linear fit
final_index = np.arange(1500,1813)
mod_fixed_unique_mag = np.delete(final_fixed_unique_mag,final_index)
mod_fixed_log_number = np.delete(final_fixed_log_number,final_index)

theoretical_slope = 0.6
theoretical_yintercept = -5

print('Theoretical linear function:')
print('Slope = %.3e' %(theoretical_slope))# Note the %.3e forces the values to be printed in scientific notation with 3 decimal places.
print('Intercept = %.3e' %(theoretical_yintercept))
print()

yerror2 = 1.0/(error_fixed_log_number)  #Error of the number counts
errors2 = 1.0/mod_fixed_log_number      #Error weightings for linear fit

fit2,cov2 = np.polyfit(mod_fixed_unique_mag,mod_fixed_log_number,1,w=1/errors2,cov=True)

measured_slope2 = fit2[0]
measured_yintercept2 = fit2[1]
sig_0_2 = np.sqrt(cov2[0,0]) #The uncertainty in the measured slope
sig_1_2 = np.sqrt(cov2[1,1]) #The uncertainty in the measured intercept

print('Fited linear function:')
print('Slope = %.3e +/- %.3e' %(measured_slope2,sig_0_2))# Note the %.3e forces the values to be printed in scientific notation with 3 decimal places.
print('Intercept = %.3e +/- %.3e' %(measured_yintercept2,sig_1_2))

plt.figure()
plt.errorbar(final_fixed_unique_mag,final_fixed_log_number,yerr=yerror2,fmt='x',capsize=2,label='Data')
plt.plot(mod_fixed_unique_mag, linear_func(mod_fixed_unique_mag,measured_slope2,measured_yintercept2),label='Linear Fit')
plt.plot(mod_fixed_unique_mag, linear_func(mod_fixed_unique_mag,theoretical_slope,theoretical_yintercept),label='Theoretical Model')

plt.xlim(7,21)
plt.ylim(-2,7)
plt.xlabel('Magnitude')
plt.ylabel('Log10(Count)')
plt.title('Number Count Plot using the Fixed Aperture Method')
#plt.xticks([200,250,300,350])
#plt.yticks([200,250,300,350])
plt.legend(loc='best')
plt.savefig('Fixed_Numbercount_Plot')
plt.show()

#%%
#Comparison of background error distributions

mean = 3420
standard_dev = 13

x_values = np.linspace(3380,3480,100000,dtype=int)

Gaussian1 = normal_distribution(mean, standard_dev, x_values)
plt.plot(x_values, Gaussian1,'rx', markersize=2,label="Global Error")

mean = 3419
standard_dev = 6

x_values = np.linspace(3380,3480,100000,dtype=int)

Gaussian2 = normal_distribution(mean, standard_dev, x_values)
plt.plot(x_values, Gaussian2,'bx', markersize=2,label="Variable Local Error")

mean = 3436
standard_dev = 10
x_values = np.linspace(3380,3480,100000,dtype=int)

Gaussian3 = normal_distribution(mean, standard_dev, x_values)
plt.plot(x_values, Gaussian3,'gx', markersize=2,label="Fixed Local Error")

plt.xlim(3380,3480)
plt.xlabel('Pixel value')
plt.ylabel('Number of Occurences')
plt.title('Comparison of the Background Error Distributions')
plt.legend()
plt.show()   

