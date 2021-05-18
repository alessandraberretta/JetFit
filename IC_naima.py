import naima
from naima.models import (ExponentialCutoffPowerLaw, PowerLaw, Synchrotron,
                          InverseCompton)
from astropy.constants import c
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import m_p


'''
def compute_B(epsB_syn, n_syn, gamma_fluid):
    return np.sqrt(32 * np.pi * np.power(m_p, 2) * epsB_syn * n_syn/u.cm**3) * gamma_fluid * c
'''
# electrons distribution
# PL = PowerLaw(1e36*u.Unit('1/eV'), 1*u.TeV, 2.5)
ECPL = ExponentialCutoffPowerLaw(1e36*u.Unit('1/eV'), 1*u.TeV, 2.1, 13*u.TeV)
SYN = Synchrotron(ECPL, B=100*u.uG)
# B = compute_B(0.001, 1, 80)
# print(B)

# synchrotron flux + photon density needed to compute the IC emission
SYN = Synchrotron(ECPL, B=100*u.uG)
Esy = np.logspace(-6, 6, 100)*u.eV
print(Esy)
Lsy = SYN.flux(Esy, distance=1*u.kpc)
print(Lsy)
R = 2
# * u.pc
phn_sy = Lsy / (4 * np.pi * R**2 * c) * 2.24


# IC (SSC) flux
IC = InverseCompton(ECPL, seed_photon_fields=[['SSC', Esy, phn_sy]])
# flux_IC = IC.flux(times, distance=0*u.Mpc)
# print(flux_IC)
spectrum_energy = np.logspace(-1, 14, 100)*u.eV
# sed_IC = IC.sed(spectrum_energy, distance=1e28*u.cm)
sed_SYN = SYN.sed(spectrum_energy, distance=1.5*u.kpc)
sed_IC = IC.sed(spectrum_energy, seed='SSC', distance=1.5*u.kpc)


# plot the figure
plt.figure(figsize=(8, 5))
plt.rc('font', family='sans')
plt.rc('mathtext', fontset='custom')
plt.loglog(spectrum_energy, sed_IC, lw=2, label='IC', c='blue')
plt.loglog(spectrum_energy, sed_SYN, lw=2, label='Sync', c='red')
plt.xlabel('Photon energy [{0}]'.format(
    spectrum_energy.unit.to_string('latex_inline')))
plt.ylabel('$E^2 dN/dE$ [{0}]'.format(
    sed_SYN.unit.to_string('latex_inline')))
plt.ylim(1e-12, 1e-5)
plt.tight_layout()
plt.legend(loc='lower left')
# plt.savefig('SSC_sed_ECPL.pdf')
plt.show()
