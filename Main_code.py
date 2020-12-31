"""
Created on Fri Nov 27 12:49:17 2020
"""
#%%
import numpy as np
from astropy.io import fits
import time

#Function that filters the data
def filt_data_function(data):
    y_len = len(data)                       #4611
    #y_len = 4611
    x_len = int(np.size(data)/y_len)        #2570
    #x_len = 2570
    #print("X length =",x_len)
    #print("Y length =",y_len)
    filt_data = data
    for i in range(x_len):
        for j in range(y_len):
        #Measurements larger than 50000 are unreliable
            if filt_data[j][i] > 50000:
                filt_data[j][i] = 0
        #Removing the edges, areas 1,2,3,4
            if i<200:
                filt_data[j][i] = 0
            if i>(x_len-200):
                filt_data[j][i] = 0
            if j<200:
                filt_data[j][i] = 0
            if j>(y_len-200):
                filt_data[j][i] = 0
        #Remove large star, area 5
            if i>1174 and i<1674:
                if j>2933 and j<3450:
                    filt_data[j][i] = 0
        #Remove large star blooming 1, area 6
            if i>1410 and i<1466:
                filt_data[j][i] = 0
        #Remove large star blooming 2, area 7
            if i>1012 and i<1707:
                if j<490:
                    filt_data[j][i] = 0
        #Remove star, area 8
            if i>726 and i<834:
                if j>3196 and j<3426:
                    filt_data[j][i] = 0
        #Remove star, area 9
            if i>923 and i<1038:
                if j>2696 and j<2839:
                    filt_data[j][i] = 0
        #Remove star, area 10
            if i>856 and i<961:
                if j>2213 and j<2367:
                    filt_data[j][i] = 0
        #Remove star, area 11
            if i>2096 and i<2182:
                if j>3703 and j<3807:
                    filt_data[j][i] = 0
        #Remove star, area 12
            if i>2174 and i<2314:
                if j>3220 and j<3345:
                    filt_data[j][i] = 0
    return filt_data

#Defines fixed aperture galaxy counter in a function
#Extend fixed aperture to include the local background calculation etc. 
def fixed_aperture(data,radius,threshold=3500):
    counting_data = data
    galaxy_count = 0
    y_len = len(counting_data)
    x_len = int(np.size(counting_data)/y_len) 
    catalogue = np.array([[0,0,0,0,0,0,0,0,0,0]])
    confirmed_wrong_counts = 0
    while np.max(counting_data) > threshold:
        max_value = np.max(counting_data)
        y_index = int(np.floor(np.argmax(counting_data)/x_len)) #Returns x-index of highest value
        x_index = int(np.argmax(counting_data)-(y_index*x_len)) #Returns y-index of highest value
        galaxy_count = galaxy_count + 1
        print('Galaxy',galaxy_count)
        print('Max value =',max_value)
        print('X Index =',x_index)
        print('Y Index =',y_index)
        print()
        #Calculating the local background
        background_lum = 0
        num_background_pixel = 0
        for j in range(int(np.floor(-2*radius)),int(np.floor(2*radius))+1):
            for k in range(int(np.floor(-2*radius)),int(np.floor(2*radius))+1):
                dist = int(np.sqrt(j**2 + k**2))
                if dist < 2 * radius:
                    if dist>radius:
                        if counting_data[y_index+j][x_index+k] != 0:
                            num_background_pixel = num_background_pixel + 1
                            background_lum = background_lum + counting_data[y_index+j][x_index+k]
        if num_background_pixel == 0:
            confirmed_wrong_counts += 1
            ave_background = 3450
        else: 
            ave_background = background_lum/num_background_pixel
        #Galaxy properties
        galaxy_lum = 0
        num_pixel = 0
        for j in range(-int(radius),int(radius)+1):
            for k in range(-int(radius),int(radius)+1):
                dist = int(np.sqrt(j**2 + k**2))
                if dist < radius:
                    add_luminosity = counting_data[y_index+j][x_index+k] - ave_background
                    if add_luminosity < 0:
                        add_luminosity=0
                    galaxy_lum = galaxy_lum + add_luminosity
                    if counting_data[y_index+j][x_index+k] != 0:
                        num_pixel = num_pixel + 1
                    counting_data[y_index+j][x_index+k] = 0
        total_lum = int(np.floor(galaxy_lum))                          
        catalogue = np.append(catalogue,[[galaxy_count,max_value,x_index,y_index,total_lum,num_pixel,radius,ave_background,2*radius,2*radius]],axis=0)
        real_catalogue = np.delete(catalogue,0,0)
    return galaxy_count, real_catalogue, confirmed_wrong_counts


hdulist = fits.open("/Users/philippsonnenschein/Desktop/Laboratory/Astro_Imaging/A1_mosaic.fits")
#hdulist = fits.open("A1_mosaic.fits")
raw_data = hdulist[0].data
filt_data = filt_data_function(raw_data)   

fixed_initial_time = time.time()
fixed_galaxy_count = fixed_aperture(filt_data,12,3500)
fixed_runtime = time.time() - fixed_initial_time

fixed_num_galaxies = fixed_galaxy_count[0]
fixed_catalogue = fixed_galaxy_count[1]
fixed_confirmed_wrong_counts = fixed_galaxy_count[2]
print("There are ",fixed_num_galaxies," Galaxies")
print("The Fixed Aperture method took",fixed_runtime,"seconds")
print()


#Variable aperture size
#Include the error of the count in the catalogue
def variable_aperture(data,threshold=3500):       #Use filtered data as input
    counting_data = data
    galaxy_count = 0
    y_len = len(counting_data)
    x_len = int(np.size(counting_data)/y_len) 
    catalogue = np.array([[0,0,0,0,0,0,0,0,0,0]])
    confirmed_wrong_counts = 0
    while np.max(counting_data) > threshold:
        max_value = np.max(counting_data)
        y_index = int(np.floor(np.argmax(counting_data)/x_len)) #Returns x-index of highest value
        x_index = int(np.argmax(counting_data)-(y_index*x_len)) #Returns y-index of highest value
        galaxy_count = galaxy_count + 1
        #finding the radius 
        rad = np.ones(4,dtype=int)
        while counting_data[y_index][x_index+rad[0]] > threshold:
            rad[0] = rad[0] + 1
        while counting_data[y_index][x_index-rad[1]] > threshold:
            rad[1] =  rad[1] + 1
        while counting_data[y_index+rad[2]][x_index] > threshold:
            rad[2] = rad[2] + 1
        while counting_data[y_index-rad[3]][x_index] > threshold:
            rad[3] = rad[3] + 1 
        radius = np.max(rad)
        #The following 2 can be used to check for the number of elliptical galaxies
        xdiameter = rad[0] + rad[1]
        ydiameter = rad[2] + rad[3]
        if radius < 4:
            radius = 4
        print('Galaxy',galaxy_count)
        print('Max value =',max_value)
        print('Radius =', radius)
        print('X diameter =', xdiameter)
        print('Y diameter =', ydiameter)
        print('X Index =',x_index)
        print('Y Index =',y_index)
        print()
        #Calculation of the local background
        #Current 2x radius is way too large, take smaller sample, maybe record what difference it makes
        background_lum = 0
        num_background_pixel = 0
        for j in range(int(np.floor(-2*radius)),int(np.floor(2*radius))+1):
            for k in range(int(np.floor(-2*radius)),int(np.floor(2*radius))+1):
                dist = int(np.sqrt(j**2 + k**2))
                if dist < 2 * radius:
                    if dist>radius:
                        if counting_data[y_index+j][x_index+k] != 0:
                            num_background_pixel = num_background_pixel + 1
                            background_lum = background_lum + counting_data[y_index+j][x_index+k]
        if num_background_pixel == 0:
            confirmed_wrong_counts += 1
            ave_background = 3450
        else: 
            ave_background = background_lum/num_background_pixel
        #Galaxy Properties
        galaxy_lum = 0
        num_pixel = 0
        for j in range(-int(radius),int(radius)+1):
            for k in range(-int(radius),int(radius)+1):
                dist = int(np.sqrt(j**2 + k**2))
                if dist < radius:
                    add_luminosity = counting_data[y_index+j][x_index+k] - ave_background
                    if add_luminosity < 0:
                        add_luminosity=0
                    galaxy_lum = galaxy_lum + add_luminosity
                    if counting_data[y_index+j][x_index+k] != 0:
                        num_pixel = num_pixel + 1
                    counting_data[y_index+j][x_index+k] = 0
        total_lum = int(np.floor(galaxy_lum))
        catalogue = np.append(catalogue,[[galaxy_count,max_value,x_index,y_index,total_lum,num_pixel,radius,ave_background,xdiameter,ydiameter]],axis=0)
        real_catalogue = np.delete(catalogue,0,0)
    return galaxy_count, real_catalogue, confirmed_wrong_counts


hdulist = fits.open("/Users/philippsonnenschein/Desktop/Laboratory/Astro_Imaging/A1_mosaic.fits")
#hdulist = fits.open("A1_mosaic.fits")
raw_data = hdulist[0].data
filt_data = filt_data_function(raw_data)  

variable_initial_time = time.time()
variable_galaxy_count = variable_aperture(filt_data,3500)
variable_runtime = time.time() - variable_initial_time

variable_num_galaxies = variable_galaxy_count[0]
variable_catalogue = variable_galaxy_count[1]
variable_confirmed_wrong_counts = variable_galaxy_count[2]
print("There are ",variable_num_galaxies," Galaxies")
print("The Variable Aperture method took",variable_runtime,"seconds")
print("Number of wrong counts:",variable_confirmed_wrong_counts)


hdulist.close("/Users/philippsonnenschein/Desktop/Laboratory/Astro_Imaging/A1_mosaic.fits")
#hdulist.close('A1_mosaic.fits') 

