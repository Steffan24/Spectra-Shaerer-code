# variables.py

from modules import np
from constants import pc, M_sun

ttt = 'ge0' #Options - ge0: no mass loss, mdt: with mass loss
imf = 'sal'
mup = '500'
low = "001"
sfh = 'is5'

n_single = 31 #SED of population (single plotting)

n_array = np.linspace(31,1030, 999, dtype='int') # SED (multiple plotting)

M_gauss = (10**6)*M_sun
d_gauss = 100*pc

save = False
n = 'single'
z = 11
R = 7000





#HARMONI spectograph resolution data (microm)

LR_IZJ_min = 0.811
LR_IZJ_max = 1.369
LR_HK_min = 1.450
LR_HK_max = 2.450

MR_IZ_min = 0.830
MR_IZ_max = 1.050
MR_J_min = 1.046
MR_J_max = 1.324
MR_H_min = 1.435
MR_H_max = 1.815
MR_K_min = 1.951
MR_K_max = 2.469
