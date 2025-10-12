### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import ascii
import latex
import os

### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 20
plt.rcParams["legend.fontsize"] = 20
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.subplot.wspace"] = 0.03
plt.rcParams["legend.frameon"] = False

### CONSTANTS USED ###
T_sun = 5772 #k
kb = 1.38*10**(-16)
c = 3*10**10
h = 6.63*10**(-27)
T_100M = 10**(4.7)
pc = 3.08*10**18
AU = 1.496*10**13
d = 1 * AU
R_sun = 6.957*10**10


### VARIABLES (unchanged as of 07/10/25)  ###
ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'


n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31,1030, 999, dtype='int') # SED (multiple plotting)


### FUNCTION FOR SINGLE LINE PLOT ####

def single_plotting(n_single):
    #read data
    
    file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
    if os.path.exists(file_loc):
        age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
        print(age)
        age_data = age['6'][0]
        print(f"AGE: {age_data}")
        age_data = (10**(age_data))*10**(-6)
        data = ascii.read(file_loc,guess = True, data_start = 2)
        wavelength = data['col1']
        total_flux = data['col3']
        total_flux = total_flux / (4*np.pi*(d**2)*(1.989*10**33))
        log_flux = np.log10(total_flux)
        log_flux = log_flux

        #PLOT 1: logf vs lambda
        plt.figure()
        plt.plot(lambda_sun, log_B+30, label = "\(Blackbody\ 1M_{\odot}\)", linestyle = '--', color = 'red', zorder = 2)
        plt.plot(wavelength, log_flux+30,c='darkblue',label = f'\({age_data} Myr\ since\ ZAMS\)', zorder = 1)
        #plt.xlim(0,7500)
        #plt.ylim(0,8)
        plt.xlim(0,6500)
        plt.xlabel("\(\lambda (\mathring{A})\)")
        plt.ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        plt.legend(bbox_to_anchor = [0.85,1], ncols = 2)
        plt.show()

        #PLOT 2: logf vs loglambda
        plt.figure()
        plt.plot(np.log10(lambda_sun), log_B+ 30, label = "\(Blackbody\ 1M_{\odot}\)", linestyle = '--', color = 'red', zorder = 2)
        plt.plot(np.log10(wavelength), log_flux + 30,c='darkblue',label = f'\({age_data} Myr\ since\ ZAMS\)', zorder = 1)
        #plt.xlim(0,7500)
        plt.ylim(0,8)
        plt.xlim(2,5.3)
        plt.xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        plt.ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        plt.legend(bbox_to_anchor = [0.85,1], ncols = 2)
        plt.show()

### FUNCTION FOR MANY PLOTS ###

def multiple_plotting(n_array):
    #PLOT 1: logf vs lambda
    fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
    ax1.plot(lambda_sun, log_B + 30, label = '\(Blackbody\ 1M_{\odot}\)', linestyle = '--', color='red', zorder = 2)
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
            total_flux = data['col3']
            total_flux = total_flux / (4*np.pi*(d**2)*(1.99*10**33))
            log_flux = np.log10(total_flux)
            ax1.plot(wavelength, log_flux + 30,c=colour[i], zorder = 1)
    all_ages_raw = []
    for i in range(len(all_ages)):
        all_ages_raw.append((10**(-6))*10**(all_ages[i]))
    ax1.set_xlim(0,7500)
    ax1.set_xlabel("\(\lambda (\mathring{A})\)")
    ax1.set_ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
    ax1.set_ylim(0,8)
    ax1.legend(bbox_to_anchor = [1.15,1.11])
    plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
    labels = [f'\({min(all_ages_raw)}\)', f'\({max(all_ages_raw)}\)']
    ax2.set_yticks([0,1030], labels=labels)
    ax2.set_ylabel('\(Time\ since\ ZAMS\ (Myr)\)')
    plt.show()

    #PLOT 2: f vs lambda
    fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
    ax1.plot(np.log10(lambda_sun), log_B + 30, label = '\(Blackbody\ 1M_{\odot}\)', linestyle = '--', color='red', zorder = 2)
    ax1.plot(np.log10(lambda_blackbody), log_blackbody + 30, label = '\(Blackbody\ 100M_{\odot}\)', linestyle = '--', color = 'orange', zorder = 2)
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
            total_flux = data['col3']
            total_flux = total_flux / (4*np.pi*(d**2)*(1.99*10**33))
            log_flux = np.log10(total_flux)
            ax1.plot(np.log10(wavelength), log_flux + 30 ,c=colour[i], zorder = 1)
    all_ages_raw = []
    for i in range(len(all_ages)):
        all_ages_raw.append((10**(-6))*10**(all_ages[i]))
    #ax1.set_xlim(0,7500)
    ax1.set_xlabel("\(\log{\lambda}\ (\mathring{A})\)")
    ax1.set_ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
    ax1.set_ylim(-2,6)
    ax1.set_xlim(2,5.2)
    ax1.legend(bbox_to_anchor = [1.15,1.15])
    plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
    labels = [f'\({min(all_ages_raw)}\)', f'\({max(all_ages_raw)}\)']
    ax2.set_yticks([0,1030], labels=labels)
    ax2.set_ylabel('\(Time\ since\ ZAMS\ (Myr)\)')
    plt.show()
        
        
### FUNCTION FOR SUN TYPE STAR ###
def sun_type_star():
    lambda_sun = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_sun_cm = lambda_sun * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_sun_cm**5))*(1/(np.exp(h*c/(lambda_sun_cm*kb*T_sun)) - 1))
    B = B * (1*10**(-8))
    B = np.pi * B * (R_sun)**2 / ((d**2)*1.99*10**33)
    log_B = np.log10(B)
    return lambda_sun, B, log_B

### FUNCTION FOR A MODEL BLACKBODY ###
def blackbody(T):
    lambda_blackbody = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_blackbody_m = lambda_blackbody * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_blackbody_m**5))*(1/(np.exp(h*c/(lambda_blackbody_m*kb*T)) - 1))
    blackbody = B * (1*10**(-8))
    blackbody = np.pi * blackbody* (R_sun*22)**2/ ((d**2)*100*1.99*10**33)
    log_blackbody = np.log10(blackbody)
    return lambda_blackbody, blackbody, log_blackbody
    
   

lambda_sun, B, log_B = sun_type_star()
lambda_blackbody, blackbody, log_blackbody = blackbody(T_100M)

single_plotting(n_single)
multiple_plotting(n_array)



