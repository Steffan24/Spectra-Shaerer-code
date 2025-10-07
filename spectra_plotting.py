### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import latex

### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 16
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True

### CONSTANTS USED ###
T_sun = 5772 #k
kb = 1.38*10**(-23)
c = 3*10**8
h = 6.63*10**(-34)

### VARIABLES (unchanged as of 07/10/25)  ###
ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'


n_single = 900 #SED of population (single plotting)

n_array = np.linspace(100,1030, 93, dtype='int') # SED (multiple plotting)


### FUNCTION FOR SINGLE LINE PLOT ####

def single_plotting(n_single):
    #read data
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
    data = ascii.read(file_loc,guess = True, data_start = 2)
    wavelength = data['col1']
    total_flux = data['col2']
    log_flux = np.log10(total_flux) #to natch schearer plot

    #PLOT 1: f vs lambda
    plt.figure()
    plt.plot(wavelength,total_flux, label = f'\({n_single}\)')
    plt.plot(lambda_sun, B, label = "\(Sun\)")
    plt.xlim(0,7500)
    plt.xlabel("\(\lambda (\mathring{A})\)")
    plt.ylabel("\(f (erg \cdot s^{-1} \cdot \mathring{A}^{-1})\)")
    plt.legend(bbox_to_anchor = [0.55,1.11], ncols = 2)
    plt.show()

### FUNCTION FOR MANY PLOTS ###

def multiple_plotting(n_array):
    plt.figure()
    plt.plot(lambda_sun, B, label = '\(Sun\)')
    for i in range(len(n_array)):
        #read dara
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
        data = ascii.read(file_loc,guess = True, data_start = 2)
        wavelength = data['col1']
        total_flux = data['col2']
        log_flux = np.log10(total_flux)
        plt.plot(wavelength, total_flux, label = f'\({n_array[i]}\)')
    plt.xlim(0,7500)
    plt.xlabel("\(\lambda (\mathring{A})\)")
    plt.ylabel("\(f (erg \cdot s^{-1} \cdot \mathring{A}^{-1})\)")
    plt.legend(bbox_to_anchor = [0.55,1.11], ncols = 2)
    plt.show()
        
        
### FUNCTION FOR SUN TYPE STAR ###
def sun_type_star():
    lambda_sun = np.linspace(0, 7500 , 1*10**5) #angstrom
    lambda_sun_m = lambda_sun * 1*10**(-10)
    B = ((2*h*c**2)/(lambda_sun_m**5))*(1/(np.exp(h*c/(lambda_sun_m*kb*T_sun)) - 1))
    log_B = np.log10(B)
    return lambda_sun, B
    
   

lambda_sun, B = sun_type_star()
single_plotting(n_single)
multiple_plotting(n_array)



