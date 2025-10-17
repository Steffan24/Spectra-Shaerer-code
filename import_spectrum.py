# IMPORT_SPECTRUM.PY (STEFFAN RHYS THOMAS L4 PROJECT - .PY FILE 1)'
# POPIII STAR MODEL VARIABLES CAN BE CHANGED DIRECTLY BELOW, EXPLANATION OF PARAMS IN README

### MODULES ###
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import ascii
import latex
import os

### VARIABLES (unchanged as of 07/10/25)  ###
ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'

n = 'multi' #single or multi

n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31,1030, 999, dtype='int') # SED (multiple plotting)

save = False #whether to save or not

T = T_100 = 10**(4.7)




### PLOTTING PARAMS ###
plt.rcParams["font.size"] = 30
plt.rcParams["legend.fontsize"] = 30
plt.rcParams["figure.frameon"] = False
plt.rcParams["figure.titlesize"] = 16
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.subplot.wspace"] = 0.03
plt.rcParams["legend.frameon"] = False
plt.rcParams["axes.linewidth"] = 1.75

### CONSTANTS USED ###
#SI UNITS
T_sun = 5772 #k'
c_m = 3*10**8 #m
T_100M = 10**(4.7) #k
M_sun_kg = 1.989*10**30
G = 6.67*10**(-11)

#non SI
kb = 1.38*10**(-16)
c = 3*10**10
h = 6.63*10**(-27)
pc = 3.08*10**18
AU = 1.496*10**13
d = 1 * AU
R_sun = 6.957*10**10
M_sun = 1.989*10**33

#stellar lifetimes
t_100 = np.log10((10**10)*(100)**(-2.5))
t_50 = np.log10((10**10)*(50)**(-2.5))
t_10 = np.log10((10**10)*(10)**(-2.5))
t_5 = np.log10((10**10)*(5)**(-2.5))

## MAIN FUNCTIONS ###

def import_data(n, save):
    #read data and determine if single or multiple files
    if n == 'single':
        file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_single}"
        if os.path.exists(file_loc):
            age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
            age_data = np.array(age['6'][0])
            data = ascii.read(file_loc,guess = True, data_start = 0)
            wavelength = np.array(data['col1'])
            flux_raw = np.array(data['col3'])
            SED_flux = flux_raw / (4*np.pi*(d**2)*M_sun)
            age_repeat = np.repeat(age_data, len(wavelength))
            SED_data = {"ages": np.array(age_repeat),
                        "wavelengths": np.array(wavelength),
                        "SED_flux": np.array(SED_flux)}
    elif n == 'multi':
        all_ages = []
        all_wavelengths = []
        all_fluxes = []
        for i in range(len(n_array)):
            file_loc = f"/home/steff/hsim/zackrisson_pop3_all/reionis_2010/pop3_{ttt}_{imf}_{mup}_{low}_{sfh}.{n_array[i]}"
            if os.path.exists(file_loc):
                age = ascii.read(file_loc,format = 'csv', comment = '$',delimiter = '\s',data_start = 12, data_end = 13, names = ['preamble1', 'preamble2','3','4','5','6',' data'])
                age_data = np.array(age['6'][0])
                data = ascii.read(file_loc,guess = True, data_start = 2)
                wavelength = np.array(data['col1'])
                all_wavelengths.append(wavelength)
                flux_raw = np.array(data['col3'])
                SED_flux = flux_raw / (4*np.pi*(d**2)*M_sun)
                age_repeat = np.repeat(age_data, len(wavelength))
                all_fluxes.append(SED_flux)
                all_ages.append(age_repeat)
        SED_data = {"ages": np.array(all_ages),
                    "wavelengths": all_wavelengths,
                    "SED_flux": all_fluxes}
    if save == True:
        np.save(f'{n}_saved_data_model_{ttt}_{imf}_{mup}_{low}_{sfh}', SED_data)
    print(SED_data)
    return SED_data

def plot(n, SED_data):
    if n == 'single':
        log_SED = np.log10(SED_data["SED_flux"])
        log_wavelength = np.log10(SED_data["wavelengths"])
        plt.figure()
        plt.plot(np.log10(lambda_sun), log_B + 30, label = "\(Blackbody\ 1M_{\odot}\)", linestyle = '--', color = 'red', zorder = 2)
        plt.plot(log_wavelength, log_SED + 30,c='darkblue',label = f'\({SED_data["ages"][0]} Myr\ since\ ZAMS\)', zorder = 1)
        #plt.xlim(0,7500)
        plt.ylim(0,8)
        plt.xlim(2,5.3)
        plt.xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        plt.ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        plt.legend(bbox_to_anchor = [0.85,1], ncols = 2)
        plt.show()
        
    elif n == 'multi':
        fig,(ax1,ax2) = plt.subplots(1,2,width_ratios=[0.95,0.05])
        ax1.plot(np.log10(lambda_sun), log_B + 30, label = '\(Blackbody\ 1M_{\odot}\)', linestyle = '--', color='red', zorder = 2)
        ax1.plot(np.log10(lambda_blackbody), log_blackbody + 30, label = '\(Blackbody\ 100M_{\odot}\)', linestyle = '--', color = 'orange', zorder = 2)
        cmap = 'viridis'
        ages = np.array([a[0] for a in SED_data["ages"]])
        #norm = Normalize(np.min(ages), np.max(ages))
        norm = Normalize(vmin=ages.min(), vmax=ages.max())
        cmap = plt.cm.viridis
        colours = cmap(norm(ages))
        bar = ScalarMappable(cmap=cmap, norm=norm)
        for i in range(len(SED_data["ages"])):
            log_wavelength = np.log10(SED_data["wavelengths"][i])
            log_SED = np.log10(SED_data["SED_flux"][i])
            age = SED_data["ages"][i][0]
            ax1.plot(log_wavelength, log_SED + 30, c=colours[i])
        ax1.set_xlabel("\(\log{\lambda}\ (\mathring{A})\)")
        ax1.set_ylabel("\(logF_{\lambda}\ 1e+30\ (erg \cdot s^{-1} \cdot \mathring{A}^{-1} \cdot cm^{-2} \cdot M_{\odot}^{-1})\)")
        ax1.set_ylim(-2,6)
        ax1.set_xlim(2,5.2)
        ax1.legend(bbox_to_anchor = [1,1])
        ax1.text(4.776,-0.38, '\(t_{100M_{\odot}}\)', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(4.709,-0.69), xycoords='data', xytext=(4.776,-0.38), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(4.4,0.2, '\(t_{50M_{\odot}}\)', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy=(4.286,0.18), xycoords='data', xytext=(4.387,0.32), textcoords = 'data', arrowprops = dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(3.994,2.68, '\(t_{10M_{\odot}}\)', bbox=dict(edgecolor='black', fc='None'))
        ax1.annotate("",xy = (3.395,0.66), xycoords='data', xytext=(3.994,2.68), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle='arc3'))
        ax1.text(4.186, 1.28, '\(t_{5M_{\odot}}\)', bbox=dict(edgecolor='black', fc = 'None'))
        ax1.annotate("",xy = (3.379,0.02), xycoords='data', xytext=(4.186,1.28), textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle='arc3'))
        plt.colorbar(bar,cax=ax2,location = 'right', orientation = 'vertical')
        labels = [f'\({min(ages)}\)','\(t_{100M_{\odot}}\)','\(t_{50M_{\odot}}\)','\(t_{10M_{\odot}}\)','\(t_{5M_{\odot}}\)', f'\({max(ages)}\)']
        ax2.set_yticks([min(ages),t_100, t_50, t_10, t_5,max(ages)], labels=labels)
        ax2.set_ylabel('\(\log{t}\ since\ ZAMS\ (yr)\)')
        plt.show()
                                    


#sun
def sun_type_star():
    lambda_sun = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_sun_cm = lambda_sun * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_sun_cm**5))*(1/(np.exp(h*c/(lambda_sun_cm*kb*T_sun)) - 1))
    B = B * (1*10**(-8))
    B = 4*np.pi * B * (R_sun)**2 /((d**2)*M_sun)
    log_B = np.log10(B)
    return lambda_sun, B, log_B


### FUNCTION FOR A MODEL BLACKBODY ###
def blackbody(T):
    lambda_blackbody = np.linspace(0, 300000 , 1*10**5) #angstrom
    lambda_blackbody_m = lambda_blackbody * 1*10**(-8)
    B = ((2*h*c**2)/(lambda_blackbody_m**5))*(1/(np.exp(h*c/(lambda_blackbody_m*kb*T)) - 1))
    blackbody = B * (1*10**(-8))
    blackbody = 4*np.pi * blackbody* (R_sun*(13.8))**2/ ((d**2)*100*M_sun)
    log_blackbody = np.log10(blackbody)
    return lambda_blackbody, blackbody, log_blackbody
    
    
## CALL FUNCTIONS

lambda_sun, B, log_B = sun_type_star()
lambda_blackbody, blackbody, log_blackbody = blackbody(T)
SED_data = import_data(n, save)
print(type(SED_data["ages"][0]), SED_data["ages"][0].shape)
print(type(SED_data["wavelengths"][0]), SED_data["wavelengths"][0].shape)
print(type(SED_data["SED_flux"][0]), SED_data["SED_flux"][0].shape)
plot(n, SED_data)
