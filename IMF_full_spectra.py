import numpy as np
import matplotlib.pyplot as plt
import plotting_params
from matplotlib.ticker import ScalarFormatter

dir = '/home/steff/hsim/HSIM/hsim/input_cubes/{model}/360ks_exposures'

imf_models = ['salpeter', 'logA', 'logE']
imf_labels = ['Salpeter', 'logL', 'logH']
adjustment = [0,100,200]


fig, ax = plt.subplots(figsize=(4, 3))
#ax = ax[0]
for i in range(len(imf_models)):
    dir = f'/home/steff/hsim/HSIM/hsim/input_cubes/{imf_models[i]}/360ks_exposures'
    wavelength = np.load(f'{dir}/M_6.0_wavelength.npy') + adjustment[i]
    flux = np.load(f'{dir}/M_6.0_flux.npy')

    ax.step(wavelength, flux, label = rf'${imf_labels[i]}$')

ax.set_xlabel('$Wavelength (\mathring{A})$')
ax.set_yscale('log')
ax.set_ylim(10**(-22), 3*10**(-17))
#ax.set_xscale('log')
ax.set_xlim(1.4*10**4, 2.2*10**4)
ax.set_ylabel('$Flux (erg \cdot s^{-1} \cdot cm^{-2} \cdot \mathring{A}^{-1})$')
ax.legend(loc = 'upper center',ncols = 3,  bbox_to_anchor =(0.5, 1.15),frameon = False,
    borderpad=0.2,
    labelspacing=0.3,
    handletextpad=0.4
)
ax.text(14750, 2*10**(-22), r'$Ly\alpha$')
ax.text(20000, 2*10**(-22), r'$HeII$')
ax.text(19300, 3*10**(-18), '$M = 10^{6.0}M_{\odot}$ \n $\Delta \lambda = 100\mathring{A}$')
plt.savefig('/home/steff/hsim/report_plots/IMF_spectra.png')
plt.show()

    
