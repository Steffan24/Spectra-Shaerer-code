import numpy as np
import matplotlib.pyplot as plt
import plotting_params

import numpy as np
import matplotlib.pyplot as plt

Mmin, Mmax = 0.5, 1000.0
M = np.logspace(np.log10(Mmin), np.log10(Mmax), 4000)

# --- Shapes in dN/dlog10M space ---
def lognormal_dndlog10M_shape(M, Mc, sigma):
    return np.exp(-(np.log(M/Mc)**2) / (2.0*sigma**2))

def salpeter_dndlog10M_shape(M, alpha=2.35):
    # dN/dM ∝ M^-alpha  ->  dN/dlog10M ∝ M^(1-alpha)
    return M**(1.0 - alpha)  # = M^-1.35

# --- Mass-normalize (same total mass formed) ---
# For dN/dlog10M: Mtot = (1/ln10) ∫ phi(M) dM
def mass_normalize(phi_dndlog10M, M, Mtot=1.0):
    integral = np.trapz(phi_dndlog10M, M)
    return phi_dndlog10M * (Mtot * np.log(10.0) / integral)

phi_sal = mass_normalize(salpeter_dndlog10M_shape(M), M, Mtot=1.0)
phi_A   = mass_normalize(lognormal_dndlog10M_shape(M, 10, 1.0), M, Mtot=1.0)
phi_E   = mass_normalize(lognormal_dndlog10M_shape(M, 60, 1.0), M, Mtot=1.0)

# --- GLOBAL SCALE to match Tumlinson y-axis placement ---
# Choose: Salpeter at 1 Msun ~ 1e-2 (from the paper panel)
target = 1e-2
scale = target / np.interp(1.0, M, phi_sal)

phi_sal *= scale
phi_A   *= scale
phi_E   *= scale

plt.figure()
plt.plot(M, phi_sal, 'k', label='$Salpeter$')
plt.plot(M, phi_A,  'r', label='$logL\ (10,1.0)$')
plt.plot(M, phi_E,  color='orange', label='$logH\ (60,1.0)$')

plt.xscale('log'); plt.yscale('log')
plt.xlim(Mmin, Mmax)
plt.ylim(1e-5, 1e-1)
plt.xlabel(r'$M\ (M_\odot)$')
plt.ylabel(r'$dN/d\log_{10}M$')
plt.legend()
plt.savefig('/home/steff/hsim/report_plots/IMF_models.png')
plt.show()
