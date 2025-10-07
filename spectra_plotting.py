### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import latex

### VARIABLES ###
ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'
n = 900 #SED of population
file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n}"

R = 7*10**8 #m
sigma = 5.67*10**(-8) # Wm^-2 K^-4
T_sun = 5772 #k
kb = 1.38*10**(-23)
c = 3*10**8
h = 6.63*10**(-34)

#print(f"FILE LOC: {file_loc}")


### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 16
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True


### IMPORT DATA ###

data = ascii.read(file_loc,guess = True, data_start = 2)


### BLACKBODY ###
Peak_lum = 4*np.pi*(R**2)*sigma*T_sun**2

lambda_sun = np.linspace(0, max(data['col1']) , 1*10**5) #angstrom
lambda_sun_m = lambda_sun * 1*10**(-10)

B = ((2*h*c**2)/(lambda_sun_m**5))*(1/(np.exp(h*c/(lambda_sun_m*kb*T_sun)) - 1))
log_B = np.log10(B)
print(lambda_sun)
print(f"B = {B}")
print(f"logB = {log_B}")

### PLOT DATA ###
wavelength = data['col1']
total_flux = data['col2']
log_flux = np.log10(total_flux)


plt.figure()
plt.plot(wavelength,total_flux, label = f'\({n}\)')
plt.plot(lambda_sun, B, label = "\(Sun\)")
plt.xlim(0,7500)
#plt.ylim(0, 35)
plt.xlabel("\(\lambda (\mathring{A})\)")
plt.ylabel("\(\log_{10}(f)\)")
#plt.title("Population III star initial model compared to blackbody")
plt.legend(bbox_to_anchor = [0.55,1.11])
plt.show()
