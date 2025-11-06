# modules.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from astropy.io import ascii
import latex
import os
from astropy.cosmology import WMAP9 as cosmo
from scipy import interpolate
import matplotlib.ticker as mticker
from astropy.table import Table
from astropy.io import fits
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

