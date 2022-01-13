import sys
import os
import numpy as np
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt
import numba
import sympy as sym

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
    R = delta*c*dt/(1+z)
    return R


def T(x):
    # Adachi \& Kasai (2012) found that T(s) can be approximated as:
    b_1 = 2.64086441
    b_2 = 0.883044401
    b_3 = 0.0531249537
    c_1 = 1.39186078
    c_2 = 0.512094674
    c_3 = 0.0394382061
    T_s = np.sqrt(x)*((2+b_1*x**3 + b_2*x**6 + b_3*x**9) /
                      (1+c_1*x**3+c_2*x**6+c_3*x**9))
    return T_s


def luminosityDistance(z):
    s = ((1.-omegaM)/omegaM)**(1./3.)
    dL = ((ckm*(1+z))/(H0*np.sqrt(s*omegaM)))*(T(s)-T(s/(1+z)))
    return dL


def doppler(gamma, theta):  # theta in rad
    beta = np.sqrt(1. - 1./(gamma*gamma))
    d = 1./(gamma*(1-beta*np.cos(theta)))
    return d


def getLogFreqArray(min, max, N):
    fre = np.logspace(min, max, num=N, endpoint=True)
    return fre


def electronDistPL(n0, gammam, gammaM, p):
    ga = np.logspace(gammam, gammaM, 200)
    ne = n0*np.power(ga, -p)
    return ne


def nu_L(b):
    k = e/(2.*np.pi*me*c)
    # print(k)
    nuL = k*b
    return nuL


def nu_s(gamma_e, B):
    nus = (4./3)*gamma_e*gamma_e*nu_L(B)
    return nus


def get_nu_i_min(gamma_e, B, gamma_min, gamma, fre):
    nus = nu_s(gamma_e, B)
    nus_min = 1.2*1e6*np.power(gamma_min, 2)*B
    # def f_nu_i_min(nu): return nus_min*(nus) / \
    # (4*np.power(gamma, 2)*(1-(h*nu)/(gamma*mec2)))
    f_nu_i_min = []
    for elm in fre:
        f_nu_i_min.append(
            nus_min*(nus)/(4*np.power(gamma, 2)*(1-(h*elm)/(gamma*mec2))))
    max_f_nu_i_min = max(f_nu_i_min)
    # max_f_nu_i_min = optimize.fmin(lambda nu: -f_nu_i_min(nu), 0)
    return max_f_nu_i_min


def get_nu_i_max(gamma_e, B, gamma_max, gamma, fre):
    nus = nu_s(gamma_e, B)
    nus_max = 1.2*1e6*np.power(gamma_max, 2)*B
    # def f_nu_i_max(nu): return nus_max*(nus)/(1-(h*nu)/(gamma*mec2))
    f_nu_i_max = []
    for elm in fre:
        f_nu_i_max.append(nus_max*(nus)/(1-(h*elm)/(gamma*mec2)))
    # max_f_nu_i_max = optimize.fmin(lambda nu: -f_nu_i_max(nu), 0)
    max_f_nu_i_max = max(f_nu_i_max)
    return max_f_nu_i_max


def flux_ssc(fre, L_ssc, d, dl):
    for elm in fre:
        F_ssc = (elm*L_ssc*np.pow(d, 4))/(4*np.pi*np.pow(dl, 2))
    return F_ssc


def luminosity_ssc(R, j_ssc):
    V_R = (4/3)*np.pi*np.pow(R, 3)
    L_ssc = (4*np.pi*V_R*j_ssc)
    return L_ssc


def emissivity_ssc(n0, gamma_min, gamma_max, p, P_ssc):
    e_dist = electronDistPL(n0, gamma_min, gamma_max, p)
    integral_j = integrate.quad(e_dist*P_ssc, gamma_min, gamma_max)
    print(integral_j)
    j_ssc = (1/4*np.pi)*integral_j
    return e_dist, j_ssc


def GAMMA_f_c(gamma, nui):
    GAMMA = (4*gamma*h*nui)/mec2
    return GAMMA


def q_f_c(gamma, nui):
    nu = 1e17  # Hz
    q = nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2)))
    return q


def int_jssc(gamma, nu_min, nu_max):
    nu = 1e17  # Hz
    integral_jssc = integrate.quad(lambda nui: (nu/(4*np.power(gamma, 2)*nui)) * (2*(nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2)))) * np.log(nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2))))+1+(nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2)))) - 2*np.power(nu/(4*np.power(gamma, 2)*nui * ( 1-(h*nu)/(gamma*mec2))), 2) + ((np.power((4*gamma*h*nui)/mec2*nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2))), 2)) * (1-nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2)))))/(2 + 2*(4*gamma*h*nui)/mec2*nu/(4*np.power(gamma, 2) * nui*(1-(h*nu)/(gamma*mec2)))))*(1/(h*nui)), nu_min, nu_max)
    return integral_jssc


def sympyintegral(A, R, nu_min, nu_max, gamma_min, gamma_max, n0, p):
    
    nui = sym.Symbol('nui')
    nu = sym.Symbol('nu')
    gamma = sym.Symbol('gamma')

    I_1 = sym.integrate((nu/4*np.power(gamma,2)*nui)*(2*nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))*(sym.log(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))*(1/(h*nui)), (nui, nu_min, nu_max))
    I_2 = sym.integrate((nu/4*np.power(gamma,2)*nui)*(1/(h*nui))*((1)+(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))-2*(sym.expand(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)**2))), (nui, nu_min, nu_max))
    I_3 = sym.integrate((nu/4*np.power(gamma,2)*nui)*(1/(h*nui))*((4*gamma*nui*h/mec2)**2)*(sym.expand((nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))**2))*(1/(2+2*(4*gamma*nui*h/mec2)*(nu/(4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))))*(1-(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))), (nui, nu_min, nu_max))
    P_ssc = A*(I_1 + I_2 + I_3)

    # P_sss_simplified = sym.simplify(A*(I_1 + I_2 + I_3))

    j_ssc_1 = (1/(4*np.pi))*sym.integrate(A*I_1*n0*(np.power(gamma, -p)),(gamma, gamma_min, gamma_max))
    j_ssc_2 = (1/(4*np.pi))*sym.integrate(A*I_2*n0*(np.power(gamma, -p)),(gamma, gamma_min, gamma_max))
    j_ssc_3 = (1/(4*np.pi))*sym.integrate(A*I_3*n0*(np.power(gamma, -p)),(gamma, gamma_min, gamma_max))
    j_ssc_tot = j_ssc_1 + j_ssc_2 + j_ssc_3

    return P_ssc, j_ssc_1, j_ssc_2, j_ssc_3, j_ssc_tot


def luminosity_ssc(R, j_ssc):
    V_R = (4/3)*np.pi*np.power(R, 3)
    L_ssc = (4*np.pi*V_R*j_ssc)
    return L_ssc


def flux_ssc(fre, L_ssc, d, dl):
    for elm in fre:
        F_ssc = (elm*L_ssc*np.power(d, 4))/(4*np.pi*np.power(dl, 2))
    return F_ssc



'''
def Jones_kernel(gamma):
    nu = 1e17
    # GAMMA = (4*gamma*h*nu_i)/mec2
    # def GAMMA_nu_i(nu_i): return (4*gamma*h*nu_i)/mec2
    def GAMMA_nu_i(nu_i): return (4*gamma*h*nu_i)/mec2
    # q = nu/(4*np.power(gamma, 2)*nu_i*(1-(h*nu)/(gamma*mec2)))
    # def q_nu_i(nu_i): return nu/(4*np.power(gamma, 2) * nu_i*(1-(h*nu)/(gamma*mec2)))

    def q_nu_i(nu_i): return nu/(4*np.power(gamma, 2)
                                 * nu_i*(1-(h*nu)/(gamma*mec2)))
    # def f_c_nu_i(nu_i): return (nu/(4*np.power(gamma, 2)*nu_i))*(2*q_nu_i(nu_i)*np.log(q_nu_i(nu_i)) + 1 + q_nu_i(nu_i) - 2 *
    # np.power(q_nu_i(nu_i), 2)+(np.power(GAMMA_nu_i(nu_i), 2)*np.power(q_nu_i(nu_i), 2)*(1-q_nu_i(nu_i)))/(2+2*GAMMA_nu_i(nu_i)*q_nu_i(nu_i)))
    def f_c_nu_i(nu_i): return (nu/(4*np.power(gamma, 2)*nu_i))*(2*q_nu_i(nu_i)*np.log(q_nu_i(nu_i)) + 1 + q_nu_i(nu_i) - 2 *
                                                                 np.power(q_nu_i(nu_i), 2)+(np.power(GAMMA_nu_i(nu_i), 2)*np.power(q_nu_i(nu_i), 2)*(1-q_nu_i(nu_i)))/(2+2*GAMMA_nu_i(nu_i)*q_nu_i(nu_i)))
    return f_c_nu_i



def power_ssc(gamma_e, B, gamma_min, gamma_max, gamma):
    u_syn = 1
    fre = getLogFreqArray(18., 24., 200)
    nu_min = get_nu_i_min(gamma_e, B, gamma_min, gamma, fre)
    nu_max = get_nu_i_max(gamma_e, B, gamma_max, gamma, fre)
    # f_c_nu_i = Jones_kernel(100)
    # for elm in f_c:
    # def f_nu_i(nu_i): return elm*(u_syn/h*nu_i)
    #def f_nu_i(nu_i): return f_c_nu_i*(u_syn/h*nu_i)
    def f_nu_i(nu_i): return lambda nu_i: f_c_nu_i*(u_syn/h*nu_i)
    integral = integrate.quad(f_nu_i, nu_min, nu_max)
    print(integral)
    P_ssc = 8*np.pi*np.power(re, 2)*c*h*integral
    return P_ssc
'''


def main():

    gamma_e = 100.
    B = 0.1
    gamma_min = 10
    gamma_max = 10000
    gamma = 100.
    dt = 1e3
    z = 1
    thetaObs = 1./gamma
    n0 = 1e4
    p = 2.5
    nu = 1e17

    d = doppler(gamma, thetaObs)
    dl = luminosityDistance(z)
    fre = getLogFreqArray(18., 24., 200)
    R = regionSize(dt, d, z)
    nu_min = get_nu_i_min(gamma_e, B, gamma_min, gamma, fre)
    nu_max = get_nu_i_max(gamma_e, B, gamma_max, gamma, fre)

    # integrale = int_jssc(gamma, nu_min, nu_max)
    A = 8*np.pi*np.power(re,2)*c*h

    P_ssc, j_ssc_1, j_ssc_2, j_ssc_3, j_ssc_tot = sympyintegral(A, R, nu_min, nu_max, gamma_min, gamma_max, n0, p)
    print('j_ssc_tot', j_ssc_tot)

    L_ssc = luminosity_ssc(R, j_ssc_tot)
    F_ssc = flux_ssc(fre, L_ssc, d, dl)

    print(F_ssc)
    print(type(F_ssc))

    

    # f_c = Jones_kernel(100)
    # P_ssc = power_ssc(gamma_e, B, gamma_min, gamma_max, gamma)
    # print(P_ssc)
    # j_ssc = emissivity_ssc(n0, gamma_min, gamma_max, p, P_ssc)
    # L_ssc = luminosity_ssc(R, j_ssc)
    # Flux_ssc = flux_ssc(fre, L_ssc, d, dl)


if __name__ == "__main__":
    main()
