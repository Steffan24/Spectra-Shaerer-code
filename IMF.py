import numpy as np
import matplotlib.pyplot as plt
import plotting_params

M = np.linspace(0.25, 10**2, 1*10**5)

#SALPETER IMF
frac_sal = M**(-2.35)

#SCALO IMF

def scalo_imf(M):
    xi = np.zeros_like(M)
    xi[M < 1] = M[M < 1] ** -1.25
    xi[(M >= 1) & (M < 10)] = M[(M >= 1) & (M < 10)] ** -2.35
    xi[M >= 10] = M[M >= 10] ** -2.7
    return xi

xi = scalo_imf(M)

## Tumlinson A

frac_A = []
M_peak = 10
for i in range(len(M)):
    frac = (1/M[i]) * np.exp(-((np.log10(M[i]) -np.log10(M_peak))**2)/(2*(1**2)))
    frac_A.append(frac)

# TUMLINSON B
frac_B = []
M_peak = 3
for i in range(len(M)):
    frac = (1/M[i]) * np.exp(-((np.log10(M[i]) -np.log10(M_peak))**2)/(2*((0.5)**2)))
    frac_B.append(frac)

# TUMLINSON E
frac_E = []
M_peak = 100
for i in range(len(M)):
    frac = (1/M[i]) * np.exp(-((np.log10(M[i]) -np.log10(M_peak))**2)/(2*(1**2)))
    frac_E.append(frac)

# LARSON
frac_larson = []
M_peak = 5
for i in range(len(M)):
    frac = (M[i]**(-1.35))*np.exp(-M_peak/M[i])
    frac_larson.append(frac)

plt.figure()
plt.plot((M), (frac_sal), c='blue', label = '\(salpeter\)')
plt.plot(M, xi, c='green', label = '\(scalo\)')
plt.plot(M, frac_A, c='red', label='\(tumlinson\ A\)')
plt.plot(M, frac_B, c='yellow', label = '\(tumlinson\ B\)')
plt.plot(M, frac_E, c='cyan', label = '\(tumlinson\ E\)')
plt.plot(M, frac_larson, c='orange', label = '\(larson\)')
plt.xscale('log', base=10)
plt.yscale('log', base=10)
plt.xlabel("\(Stellar\ mass\ (M_{\odot})\)")
plt.ylabel(r"$\frac{dN}{dM}$")
plt.legend()
plt.show()
