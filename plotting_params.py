# plotting+params.py

from modules import plt
from modules import  np
import scienceplots

# plt.rcParams["font.size"] = 30
# plt.rcParams["legend.fontsize"] = 30
# plt.rcParams["figure.frameon"] = False
# plt.rcParams["figure.titlesize"] = 16
# plt.rcParams["text.usetex"] = True
# plt.rcParams["figure.subplot.wspace"] = 0.03
# plt.rcParams["legend.frameon"] = False
# plt.rcParams["axes.linewidth"] = 1.75
# plt.rcParams["figure.subplot.hspace"] = 0
# plt.rcParams['figure.figsize'] = [10,6]
# plt.rcParams['lines.marker'] = "."  #None
# plt.rcParams['lines.markersize'] = 6   #6.0
# plt.rcParams['lines.linewidth'] = 3 # 1.5

plt.style.use(['science', 'nature'])

plot_params = {
    "figure.dpi": "200",
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "mathtext.default": "regular",
    "axes.labelsize": 12,
    "axes.linewidth": 1.5,
    "axes.titlesize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.title_fontsize": 12,
    "legend.fontsize": 12,
    "xtick.major.size": 3.5,
    "xtick.major.width": 1.5,
    "xtick.minor.size": 2.5,
    "xtick.minor.width": 1.5,
    "ytick.major.size": 3.5,
    "ytick.major.width": 1.5,
    "ytick.minor.size": 2.5,
    "ytick.minor.width": 1.5,
    "font.size": 12
}

plt.rcParams.update(plot_params)
