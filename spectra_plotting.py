### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import ascii
import latex
import os

### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 16
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.subplot.wspace"] = 0.03

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


n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31,1030, 93, dtype='int') # SED (multiple plotting)


### FUNCTION FOR SINGLE LINE PLOT ####

def single_plotting(n_single):
    #read data
    
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
    if os.path.exists(file_loc):
        age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
        print(age)
        age_data = age['6'][0]
        print(f"AGE: {age_data}")
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

        #PLOT 2: logf vs lambda
        plt.figure()
        plt.plot(wavelength, np.log(total_flux), label = f'\({n_single}\)')
        plt.plot(lambda_sun, np.log(B), label = '\(Sun\)')
        plt.xlim(0,7500)
        plt.xlabel('\(\lambda (\mathring{A})\)')
        plt.ylabel('\(\log_{10}{f}\)')
        plt.legend(bbox_to_anchor = [0.55,1.11], ncols = 2)
        plt.show()

### FUNCTION FOR MANY PLOTS ###

def multiple_plotting(n_array):
    #PLOT 1: f vs lambda
    fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
    log_B = np.log(B) / (1.99*10**30)
    ax1.plot(lambda_sun, log_B*(10**30), label = '\(Sun\)', linestyle = '--', color='red')
    colour = plt.cm.viridis(np.linspace(0,1,len(n_array)))
    cmap = 'viridis'
    norm = Normalize(0,max(n_array))
    bar = ScalarMappable(cmap=cmap, norm=norm)
    all_ages = []
    for i in range(len(n_array)):
        #read dara
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
        if os.path.exists(file_loc):
            age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
            print(age)
            age_data = age['6'][0]
            all_ages.append(age_data)
            data = ascii.read(file_loc,guess = True, data_start = 2)
            wavelength = data['col1']
            total_flux = data['col2']
            log_flux = np.log10(total_flux)
            log_flux = log_flux / (1.99*10**30)
            ax1.plot(wavelength, log_flux*(10**30),c=colour[i])
    all_ages_raw = []
    for i in range(len(all_ages)):
        all_ages_raw.append((10**(-6))*10**(all_ages[i]))
    ax1.set_xlim(0,7500)
    ax1.set_xlabel("\(\lambda (\mathring{A})\)")
    ax1.set_ylabel("\(logF_{\lambda} 1e+30 (erg \cdot s^{-1} \cdot \mathring{A}^{-1} M_{o}^{-1})\)")
    ax1.set_ylim(0,35)
    ax1.legend(bbox_to_anchor = [0.7,1.11])
    plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
    labels = [min(all_ages_raw), max(all_ages_raw)]
    ax2.set_yticks([0,1030], labels=labels)
    ax2.set_ylabel('\(Stellar age (Myr)\)')
    plt.show()

    #PLOT 2: logf vs lambda
    plt.figure()
    plt.plot(lambda_sun, np.log10(B), label = '\(Sun\)')
    for i in range(len(n_array)):
        #read dara
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
        if os.path.exists(file_loc):
            data = ascii.read(file_loc,guess = True, data_start = 2)
            wavelength = data['col1']
            total_flux = data['col2']
            log_flux = np.log(total_flux)
            plt.plot(wavelength, log_flux, label = f'\({n_array[i]}\)')
    plt.xlim(0,7500)
    plt.ylim(0,500)
    plt.xlabel('\(\lambda (\mathring{A})\)')
    plt.ylabel('\(\log_{10}{f}\)')
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



