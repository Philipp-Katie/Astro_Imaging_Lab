"""
Created on Wed Dec 23 13:29:19 2020
"""
#%%
#Split up data into 8 segments

x1 = 1285
x2 = 2569
y1 = 1152
y2 = 2305
y3 = 3458
y4 = 4610

filt_data = filt_data_function(raw_data)  

data1 = filt_data[:y1,:x1]
data2 = filt_data[:y1,x1:x2]
data3 = filt_data[y1:y2,:x1]
data4 = filt_data[y1:y2,x1:x2]
data5 = filt_data[y2:y3,:x1]
data6 = filt_data[y2:y3,x1:x2]
data7 = filt_data[y3:y4,:x1]
data8 = filt_data[y3:y4,x1:x2]
#%% 
#Testing the counting galaxies algorithm
#Should return 0 galaxies

test_zeros=np.zeros((2570,4611))            #make array of just zeros

zero_galaxies = variable_aperture(test_zeros,3500)

fake_num_galaxies = zero_galaxies[0]
fake_catalogue = zero_galaxies[1]
fake_confirmed_wrong_counts = zero_galaxies[2]
print(fake_num_galaxies)

#%%
#r#Testing the counting galaxies algorithm using a random noise distribution

random_noise = 3450 + (60*np.random.rand(4611,2570))       #make array of random noise

filt_random_noise = filt_data(random_noise)

random_galaxies = variable_aperture(filt_random_noise,3500)

fake_num_galaxies = zero_galaxies[0]
fake_catalogue = zero_galaxies[1]
fake_confirmed_wrong_counts = zero_galaxies[2]
print(fake_num_galaxies)

#%%
#Example of the error caused by the Fixed Aperture method

#make some fake gaussian shaped data
def Gamma2sigma(Gamma):
    '''Function to convert FWHM (Gamma) to standard deviation (sigma)'''
    return Gamma * np.sqrt(2) / ( np.sqrt(2 * np.log(2)) * 2 )

x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

mu_x1 = 1255
mu_y1 = 2305
sigmax1 = Gamma2sigma(12)
z1 =  np.exp(-((x-mu_x1)**2+(y-mu_y1)**2)/(2*sigmax1**2))
z1_scaled = z1*40000  

z_scaled = z1_scaled 


plt.contourf(x, y, z_scaled, cmap='Greys')
plt.colorbar()
plt.axis('equal')
plt.show()


test_gaussian = fixed_aperture(z_scaled,12,3500)
num_fake_galaxies = test_gaussian[0]
fake_catalogue = test_gaussian[1]
print(num_fake_galaxies)

size = len(fake_catalogue)
fake_x_coordinates = np.zeros(size)
fake_y_coordinates = np.zeros(size)
fake_radius = np.zeros(size)
for i in range(0,size):
    fake_x_coordinates[i] = fake_catalogue[i][2]
    fake_y_coordinates[i] = fake_catalogue[i][3]
    fake_radius[i] = fake_catalogue[i][6]
#plt.figure()
plt.axis('equal')
plt.xlim(0,2600)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Detected Galaxies')
for i in range(size):
    plt.scatter(fake_x_coordinates[i],fake_y_coordinates[i],c='b',marker='o',s=int(15*fake_radius[i]**2))
plt.show()

#%%
#Example of the error caused by the Variable Aperture method

#make some fake gaussian shaped data
x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

mu_x1 = 1255
mu_y1 = 2305
sigmax1 = Gamma2sigma(12)
z1 =  np.exp(-((x-mu_x1)**2+(y-mu_y1)**2)/(2*sigmax1**2))
z1_scaled = z1*20000 

mu_x2 = 1245
mu_y2 = 2315
sigmax2 = Gamma2sigma(12)
z2 =  np.exp(-((x-mu_x2)**2+(y-mu_y2)**2)/(2*sigmax2**2))
z2_scaled = z2*20000


z_scaled = z1_scaled + z2_scaled


plt.contourf(x, y, z_scaled, cmap='Greys')
plt.colorbar()
plt.axis('equal')
plt.show()


test_gaussian = variable_aperture(z_scaled,3500)
num_fake_galaxies = test_gaussian[0]
fake_catalogue = test_gaussian[1]
print(num_fake_galaxies)

size = len(fake_catalogue)
fake_x_coordinates = np.zeros(size)
fake_y_coordinates = np.zeros(size)
fake_radius = np.zeros(size)
for i in range(0,size):
    fake_x_coordinates[i] = fake_catalogue[i][2]
    fake_y_coordinates[i] = fake_catalogue[i][3]
    fake_radius[i] = fake_catalogue[i][6]
#plt.figure()
plt.axis('equal')
plt.xlim(0,2600)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Detected Galaxies')
for i in range(size):
    plt.scatter(fake_x_coordinates[i],fake_y_coordinates[i],c='b',marker='o',s=int(15*fake_radius[i]**2))
plt.show()

#%%
#Example of the error caused by the Variable Aperture method

#make some fake gaussian shaped data
x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

mu_x1 = 1255
mu_y1 = 2305
sigmax1 = Gamma2sigma(12)
z1 =  np.exp(-((x-mu_x1)**2+(y-mu_y1)**2)/(2*sigmax1**2))
z1_scaled = z1*10000 

mu_x2 = 1240
mu_y2 = 2320
sigmax2 = Gamma2sigma(12)
z2 =  np.exp(-((x-mu_x2)**2+(y-mu_y2)**2)/(2*sigmax2**2))
z2_scaled = z2*10000


z_scaled = z1_scaled + z2_scaled


plt.contourf(x, y, z_scaled, cmap='Greys')
plt.colorbar()
plt.axis('equal')
plt.show()


test_gaussian = variable_aperture(z_scaled,3500)
num_fake_galaxies = test_gaussian[0]
fake_catalogue = test_gaussian[1]
print(num_fake_galaxies)

size = len(fake_catalogue)
fake_x_coordinates = np.zeros(size)
fake_y_coordinates = np.zeros(size)
fake_radius = np.zeros(size)
for i in range(0,size):
    fake_x_coordinates[i] = fake_catalogue[i][2]
    fake_y_coordinates[i] = fake_catalogue[i][3]
    fake_radius[i] = fake_catalogue[i][6]
#plt.figure()
plt.axis('equal')
plt.xlim(0,2600)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Locations of Detected Galaxies')
for i in range(size):
    plt.scatter(fake_x_coordinates[i],fake_y_coordinates[i],c='b',marker='o',s=int(15*fake_radius[i]**2))
plt.show()


#%%
#Creating more complicated fake images

x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

mu_x1 = 1255
mu_y1 = 2305
sigmax1 = Gamma2sigma(12)
z1 =  np.exp(-((x-mu_x1)**2+(y-mu_y1)**2)/(2*sigmax1**2))
z1_scaled = z1*40000  

mu_x2 = 1220
mu_y2 = 2305
sigmax2 = Gamma2sigma(12)
z2 =  np.exp(-((x-mu_x2)**2+(y-mu_y2)**2)/(2*sigmax2**2))
z2_scaled = z2*20000

mu_x3 = 1190
mu_y3 = 2305
sigmax3 = Gamma2sigma(12)
z3 =  np.exp(-((x-mu_x3)**2+(y-mu_y3)**2)/(2*sigmax3**2))
z3_scaled = z3*10000

mu_x4 = 1165
mu_y4 = 2305
sigmax4 = Gamma2sigma(12)
z4 =  np.exp(-((x-mu_x4)**2+(y-mu_y4)**2)/(2*sigmax4**2))
z4_scaled = z4*5000


z_scaled = z1_scaled + z2_scaled + z3_scaled + z4_scaled


plt.contourf(x, y, z_scaled, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()


test_gaussian = variable_aperture(z_scaled,3500)
num_fake_galaxies = test_gaussian[0]
print(num_fake_galaxies)

#%%
#Creating more sophisticated random images and comparing the 2 methods
#This uses circular galaxies

#Returns a random number betweeen x0 and x1, uniformly distributed 

import matplotlib.pyplot as plt

def random_number(x0,x1):
    return int(x0 + ((x1-x0) * np.random.rand()))

def create_circular_data(num_galaxies):
    number_of_galaxies_created = num_galaxies
    #number_of_galaxies_created = random_number(10,30)
    x = np.linspace(0,2570,2570)
    y = np.linspace(0, 4611, 4611)
    x, y = np.meshgrid(x, y)
    size = number_of_galaxies_created
    mu_x = random_number(200,2370)
    mu_y = random_number(200,4411)
    sigma = random_number(4,15)
    z_scale = random_number(3500,10000)
    layer1 = np.array(z_scale*np.exp(-((x-mu_x)**2+(y-mu_y)**2)/(2*sigma**2)))
    real_layer = np.reshape(layer1,(1,4611,2570))
    z_final = real_layer
    mu_x = np.zeros(size)
    mu_y = np.zeros(size)
    sigma = np.zeros(size)
    z_scale = np.zeros(size)
    z = real_layer
    for i in range(size):
        print("Galaxies created: ",i)
        mu_x[i] = random_number(200,2370)
        mu_y[i] = random_number(200,4411)
        sigma[i] = random_number(4,15)
        z_scale[i] = random_number(3500,10000)
        layer = np.array(z_scale[i]*np.exp(-((x-mu_x[i])**2+(y-mu_y[i])**2)/(2*sigma[i]**2)))
        real_layer = np.reshape(layer,(1,4611,2570))
        z = np.concatenate((z,real_layer), axis = 0)
    for i in range(size):
        z_final = z_final + z[i]
    z_final = np.reshape(z_final,(4611,2570))
    return z_final

#Set the number of galaxies to be created
num_galaxy = 20

x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

data1 = create_circular_data(num_galaxy)

plt.figure()
plt.contourf(x, y, data1, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()

fake1_count = fixed_aperture(data1,12,3500)
fake1_num_galaxies_counted = fake1_count[0]
fake1_catalogue = fake1_count[1]
fake1_confirmed_wrong_counts = fake1_count[2]
difference1 = abs(num_galaxy - fake1_num_galaxies_counted)
error1 = (difference1/num_galaxy)*100
print("Number of galaxies created = ",num_galaxy)
print()
print("Counted number of galaxies = ",fake1_num_galaxies_counted)
print()
print("Difference in galaxies = ", difference1)
print("The error of the fixed aperture method is",error1,"%")

data2 = create_circular_data(num_galaxy)

plt.figure()
plt.contourf(x, y, data2, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()

fake2_count = variable_aperture(data2,3500)
fake2_num_galaxies_counted = fake2_count[0]
fake2_catalogue = fake2_count[1]
fake2_confirmed_wrong_counts = fake2_count[2]
difference2 = abs(num_galaxy - fake2_num_galaxies_counted)
error2 = (difference2/num_galaxy)*100
print("Number of galaxies created = ",num_galaxy)
print()
print("Counted number of galaxies = ",fake2_num_galaxies_counted)
print()
print("Difference in galaxies = ", difference2)
print("The error of the variable aperture method is",error2,"%")
print()

if error1 > error2:
    print("The variable aperture method was",error1-error2,"% more accurate")
else:
    print("The fixed aperture method was",error2-error1,"% more accurate")
#%%
#Plotting the created fake data after the counting has been done 
#Demonstrates how the algorithms work and displays errors
plt.figure()
plt.contourf(x, y, data1, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()
#%%
plt.figure()
plt.contourf(x, y, data2, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()
#%%
#Code used to create elliptical galaxies:
    
x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

mu_x_elliptical = 1255
mu_y_elliptical = 2305
sigmax_elliptical = Gamma2sigma(5)
sigmay_elliptical = Gamma2sigma(50)
z_elliptical =  np.exp(-((((x-mu_x_elliptical)*sigmay_elliptical)**2 + ((y-mu_y_elliptical)*sigmax_elliptical)**2)/(2*sigmax_elliptical*sigmay_elliptical)**2))

z_scaled_elliptical = z_elliptical*40000  

plt.contourf(x, y, z_scaled_elliptical, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()

#%%
#Creating more sophisticated random images and comparing the 2 methods
#This uses eliptical galaxies

def create_elliptical_data(num_galaxies):
    number_of_galaxies_created = num_galaxies
    #number_of_galaxies_created = random_number(10,30)
    x = np.linspace(0,2570,2570)
    y = np.linspace(0, 4611, 4611)
    x, y = np.meshgrid(x, y)
    size = number_of_galaxies_created
    mu_x = random_number(200,2370)
    mu_y = random_number(200,4411)
    sigmax = random_number(4,15)
    sigmay = random_number(4,15)
    z_scale = random_number(3500,10000)
    layer1 = np.array(z_scale*np.exp(-((((x-mu_x)*sigmay)**2 + ((y-mu_y)*sigmax)**2)/(2*sigmax*sigmay)**2)))
    real_layer = np.reshape(layer1,(1,4611,2570))
    z_final = real_layer
    print("Galaxies created: ",1)
    mu_x = np.zeros(size)
    mu_y = np.zeros(size)
    sigmax = np.zeros(size)
    sigmay = np.zeros(size)
    z_scale = np.zeros(size)
    z = real_layer
    for i in range(1,size):
        print("Galaxies created: ",i+1)
        mu_x[i] = random_number(200,2370)
        mu_y[i] = random_number(200,4411)
        sigmax[i] = random_number(4,15)
        sigmay[i] = random_number(4,15)
        z_scale[i] = random_number(3500,10000)
        layer = np.array(z_scale[i]*np.exp(-((((x-mu_x[i])*sigmay[i])**2 + ((y-mu_y[i])*sigmax[i])**2)/(2*sigmax[i]*sigmay[i])**2)))
        real_layer = np.reshape(layer,(1,4611,2570))
        z = np.concatenate((z,real_layer), axis = 0)
    for i in range(size):
        z_final = z_final + z[i]
    z_final = np.reshape(z_final,(4611,2570))
    return z_final

#Set the number of galaxies to be created
num_galaxy = 20

data3 = create_elliptical_data(num_galaxy)
fake3_count = fixed_aperture(data3,12,3500)
fake3_num_galaxies_counted = fake3_count[0]
fake3_catalogue = fake3_count[1]
fake3_confirmed_wrong_counts = fake3_count[2]
difference3 = abs(num_galaxy - fake3_num_galaxies_counted)
error3 = (difference3/num_galaxy)*100
print("Number of galaxies created = ",num_galaxy)
print()
print("Counted number of galaxies = ",fake3_num_galaxies_counted)
print()
print("Difference in galaxies = ", difference3)
print("The error of the fixed aperture method is",error3,"%")

data4 = create_elliptical_data(num_galaxy)
fake4_count = variable_aperture(data4,3500)
fake4_num_galaxies_counted = fake4_count[0]
fake4_catalogue = fake4_count[1]
fake4_confirmed_wrong_counts = fake4_count[2]
difference4 = abs(num_galaxy - fake4_num_galaxies_counted)
error4 = (difference4/num_galaxy)*100
print("Number of galaxies created = ",num_galaxy)
print()
print("Counted number of galaxies = ",fake4_num_galaxies_counted)
print()
print("Difference in galaxies = ", difference4)
print("The error of the variable aperture method is",error4,"%")
print()

if error3 > error4:
    print("The variable aperture method was",error3-error4,"% more accurate")
else:
    print("The fixed aperture method was",error4-error3,"% more accurate")
#%%
#Plotting the created fake data after the counting has been done 
#Demonstrates how the algorithms work and displays errors

x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

plt.contourf(x, y, data3, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()
#%%
x = np.linspace(0,2570,2570)
y = np.linspace(0, 4611, 4611)
x, y = np.meshgrid(x, y)

plt.contourf(x, y, data4, cmap='Blues')
plt.colorbar()
plt.axis('equal')
plt.show()
