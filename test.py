import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def gaussian(x,A):
    std = 3.244
    B = 23207.3
    gaussian = (A/(np.sqrt(2*np.pi)*std))*np.exp(-0.5*((x-B)**2)/std**2)
    return gaussian


def calcSNR_chi2(xx, yy):
    popt, pcov = curve_fit(gaussian, xx, yy, p0=[1])


def calcSNR_peak(xx, yy):
    popt, pcov = curve_fit(gaussian, xx, yy, p0=[1])
    A = popt[0]
    peak = A / np.sqrt(2*np.pi) / 3.24
    std = np.nanstd(yy[xx > 23300])
    return peak / std * np.sqrt(2.355 * 3.244 / 3.3)


xx = np.linspace(23000, 23500, 150)
yy = np.random.normal(0, .1, 150)
yy += gaussian(xx, 3)

print(calcSNR_peak(xx, yy))

plt.step(xx, yy)
plt.xlim(23160, 23250)
#plt.show()
