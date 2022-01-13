import sys
import os
from math import *
import numpy as np
from scipy import integrate
import scipy.special as sc
import matplotlib.pyplot as plt

mec2 = 8.187e-7  # ergs
mec2eV = 5.11e5  # eV
mpc2 = 1.5032e-3  # ergs
eV2Hz = 2.418e14
eV2erg = 1.602e12
kB = 1.3807e-16  # [erg/K]
h = 6.6261e-27  # erg*sec
me = 9.1094e-28  # g
mp = 1.6726e-24  # g
G = 6.6726e-8  # dyne cm^2/g^2
Msun = 1.989e33  # g
Lsun = 3.826e33  # erg/s
Rsun = 6.960e10  # cm
pc = 3.085678e18  # cm
e = 4.8032e-10  # statcoulonb
re = 2.8179e-13  # cm
sigmaT = 6.6525e-25  # cm^2
sigmaSB = 5.6705e-5  # erg/(cm^2 s K^4 )
Bcr = 4.414e13  # G
c = 2.99792e10  # cm/s
ckm = 299792  # km/s
H0 = 67.4  # km s^-1 Mpc^-1 (Planck Collaboration 2018)
omegaM = 0.315  # (Planck Collaboration 2018)


def regionSize(dt, delta, z):
    # given doppler variation timescale (sec) and redshift
    R = delta*c*dt/(1+z)
    return R


def doppler(gamma, theta):  # theta in rad
    beta = sqrt(1. - 1./(gamma*gamma))
    d = 1./(gamma*(1-beta*cos(theta)))
    return d


def bulkGamma(delta):
    g = 1./(delta*(1-beta*cos(theta)))
    return g


def regionLocation(gamma, dt, z):
    z0 = 2.*gamma*gamma*c*dt/(1+z)
    return z0


def T(x):
    # Adachi \& Kasai (2012) found that T(s) can be approximated as:
    b_1 = 2.64086441
    b_2 = 0.883044401
    b_3 = 0.0531249537
    c_1 = 1.39186078
    c_2 = 0.512094674
    c_3 = 0.0394382061
    T_s = sqrt(x)*((2+b_1*x**3 + b_2*x**6 + b_3*x**9) /
                   (1+c_1*x**3+c_2*x**6+c_3*x**9))
    return T_s


def luminosityDistance(z):
    s = ((1.-omegaM)/omegaM)**(1./3.)
    dL = ((ckm*(1+z))/(H0*sqrt(s*omegaM)))*(T(s)-T(s/(1+z)))
    return dL


def nu_L(b):
    k = e/(2.*pi*me*c)
    # print(k)
    nuL = k*b
    return nuL


def nu_s(gamma_e, B):
    nus = (4./3)*gamma_e*gamma_e*nu_L(B)
    return nus


def nu_c(gamma_e, B):
    nus = (3./2)*gamma_e*gamma_e*nu_L(B)
    return nus


def uB(b):  # magnetic energy density
    ub = b**2/(8*pi)
    return ub


def getLogGammaArray(min, max, N):
    ar = np.logspace(min, max, num=N, endpoint=True)
    return ar


def getLinGammaArray(min, max, N):
    ar = np.linspace(min, max, num=N, endpoint=True)
    return ar


def getLogFreqArray(min, max, N):
    fre = np.logspace(min, max, num=N, endpoint=True)
    return fre


def getLinFreqArray(min, max, N):
    fre = np.linspace(min, max, num=N, endpoint=True)
    return fre


def syncEmissivityDeltaPL(freq, B, n0, ind):
    k = 3*sigmaT*c*uB(B)*n0/8.
    k1 = pow(nu_L(B), (ind-3.)/2.)
    j = k*k1*np.power(freq, -(ind-1)/2.)
    return j


def bessel(xs):
    o = 5./3.
    def bes(x): return sc.kv(o, x)
    r = integrate.quad(bes, xs, np.inf, limit=100)
    #print (xs)
    return r[0]*xs


def plotBessel():
    x1 = np.linspace(0.01, 1, 9, endpoint=True)
    o = 5./3
    y = []
    for xs in x1:
        def bes(x): return sc.kv(o, x)
        r = integrate.quad(bes, xs, 100., limit=300)
        y.append(xs*r[0])
        print(xs, y)
    print
    plt.plot(np.log10(x1), np.log10(y))
    plt.ylim(-2., 1.)
    # plt.legend()
    plt.show()


def electronDistPL(n0, gammam, gammaM, p):
    ga = getLogGammaArray(gammam, gammaM, 200)
    ne = n0*np.power(ga, -p)
    return ga, ne


def singleElecSynchroPower(nu, gamma, B):
    nuL = nu_L(B)
    n1 = 2.*pi*sqrt(3.)*e*e*nuL/c
    nus = nu_c(gamma, B)
    x0 = nu/nus
    y0 = bessel(x0)
    P = y0
    #print ("1--------------->",gamma,x0)
    return P


def syncEmissivityKernPL(gamma, nu, p, B):
    ne = pow(gamma, -p)
    k1 = ne*singleElecSynchroPower(nu, gamma, B)
    # print(ne,k1,nu,p,B,gamma)
    return k1


def syncEmissivity(freq, gammam, gammaM, p, n0, B, R):
    nuL = nu_L(B)
    f1 = []
    y = []
    assorb = []
    n1 = 2.*pi*sqrt(3.)*e*e*nuL/c
    k0 = (p+2)/(8*pi*me)
    ar = getLogGammaArray(log10(gammam), log10(gammaM), 100)
    for f in freq:
        y1 = []
        yy = []
        for x in ar:
            asso = singleElecSynchroPower(f, x, B)*pow(x, -(p+1))
            js = syncEmissivityKernPL(x, f, p, B)
            y1.append(js)
            yy.append(asso)
            # print("------>",x,yy)
        r = integrate.simps(y1)
        al = integrate.simps(yy)
        #print (f,r)
        if r > 1e-50:
            I = n1*r*n0  # (r[0]/alpha)*(1.-exp(-tau))
            as1 = n1*al*k0*n0/pow(f, 2.)
            tau = as1*R
            if(tau > 0.1):
                assorb.append(f*I*R)
                I = I*R*(1.-exp(-tau))/tau
                y.append(f*I)
                # y1.append(alpha)
                f1.append(f)
            else:
                y.append(f*I*R)
                assorb.append(f*I*R)
                # y1.append(alpha)
                f1.append(f)
    print(f1)
    plt.plot(np.log10(f1), np.log10(y))
    # plt.plot(np.log10(f1), np.log10(assorb))
    # plt.plot(np.log10(f1), np.log10(y1))
    # plt.legend()
    plt.ylim(-10, 10)
    plt.xlim(7., 25.)
    plt.show()


def syncAbsorption(freq, gammam, gammaM, p, n0, B):
    k0 = (p+2)*n0/(8.*pi*me)
    y = []
    f1 = []

    for f in freq:
        def kee(x): return singleElecSynchroPower(f, x, B)*pow(x, -(p+1))
        r = integrate.quad(kee, gammam, gammaM, limit=100)
        if r[0] > 1e-40:
            y.append(k0*r[0]/pow(f, 2.))
            f1.append(f)
    # print(y)
    plt.plot(np.log10(f1), np.log10(f1*y))
    plt.xlabel("Log(nu)")
    plt.ylabel("Log(nu*F(nu))")
    # plt.legend()
    plt.ylim(-40, 0)
    plt.show()


'''
def Jones_kernel(gamma, nu_i, fre, q, GAMMA):
    for nu in fre:
        GAMMA = (4*gamma*h*nu)/mec2
        q = nu/(4*np.power(gamma, 2)*nu*(1-(h*nu)/(gamma*mec2)))
        f_c = fre/(4*np.power(gamma, 2)*nu_i)*(2*q*np.log(q) + 1 + q - 2 *
                                               np.power(q, 2)+(np.power(GAMMA, 2)*np.power(q, 2)*(1-q))/(2+2*GAMMA*q))


def power_ssc():

    integral = integrate.quad()
    P_ssc = 8*np.pi*np.power(re, 2)*c*h*integral


def emissivity_ssc(n0, gamma_min, gamma_max, p, P_ssc):
    e_dist = electronDistPL(n0, gamma_min, gamma_max, p)
    integral_j = integrate.quad(e_dist*P_ssc, gamma_min, gamma_max)
    j_ssc = (1/4*np.pi)*integral_j


def luminosity_ssc(R, j_ssc):
    V_R = (4/3)*np.pi*np.pow(R, 3)
    L_ssc = (4*np.pi*V_R*j_ssc)


def flux_ssc(fre, L_ssc, d, dl):
    for elm in fre:
        F_ssc = (elm*L_ssc*np.pow(d, 4))/(4*np.pi*np.pow(dl, 2))
    return F_ssc
'''

if __name__ == "__main__":
    gammaL = 100.  # bulk Lorentz
    dt = 1e3  # sec
    z = 1  # redshift
    thetaObs = 1./gammaL
    gamma_e = 100.  # lorentz factor of an electron
    B = 0.1  # gauss magnetic field
    gamma_min = 10.
    gamma_max = 100000.
    n0 = 10e4
    p = 2.5
    d = doppler(gammaL, thetaObs)
    # print(d)
    R = regionSize(dt, d, z)
    # print(R, "cm", R/pc, "pc")
    dz = regionLocation(gammaL, dt, z)
    # print(dz, "cm", dz/pc, "pc")
    dl = luminosityDistance(z)
    # print(dl, "Mpc", dl*pc, "cm")
    nus = nu_s(10, B)
    # print(nus, "Hz")
    fre = getLogFreqArray(5., 18., 400)
    # print(fre)
    # plotBessel()
    syncEmissivity(fre, gamma_min, gamma_max, p, n0, B, R)
