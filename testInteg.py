import sys,os
from math import *
import numpy as np
from scipy import integrate
import scipy.special as sc

mec2 = 8.187e-7 #ergs
mec2eV = 5.11e5 #eV
mpc2 = 1.5032e-3 #ergs
eV2Hz = 2.418e14
eV2erg = 1.602e12
kB = 1.3807e-16 #[erg/K]
h = 6.6261e-27 #erg*sec
me = 9.1094e-28 #g
mp = 1.6726e-24 #g
G = 6.6726e-8 #dyne cm^2/g^2
Msun = 1.989e33 #g
Lsun = 3.826e33 #erg/s
Rsun = 6.960e10 #cm
pc = 3.085678e18 #cm
e = 4.8032e-10 #statcoulonb
re = 2.8179e-13 #cm
sigmaT = 6.6525e-25 #cm^2
sigmaSB = 5.6705e-5 # erg/(cm^2 s K^4 )
Bcr = 4.414e13 # G
c=2.99792e10 #cm/s
ckm=299792 #km/s

print(0.29*3*e/(4*pi*me*c))
