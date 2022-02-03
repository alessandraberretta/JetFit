from cmath import isnan
import os
import sys
import numpy as np
from scipy import integrate
from scipy import optimize
import scipy.special as sc
import matplotlib.pyplot as plt
import mpmath
import sympy as sym
from sympy import besselk
from sympy import oo
from sympy import sqrt
from sympy.functions import exp
from tqdm import tqdm

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

def getLogGammaArray(min, max, N):
    ar = np.logspace(min, max, num=N, endpoint=True)
    return ar

def electronDistPL(n0, gammam, gammaM, p):
    ga = np.logspace(gammam, gammaM, 200)
    ne = n0*np.power(ga, -p)
    return ne

def nu_L(b):
    k = e/(2.*np.pi*me*c)
    nuL = k*b
    return nuL

def bessel(xs):
    o = 5./3.
    def bes(x): return sc.kv(o, x)
    r = integrate.quad(bes, xs, np.inf, limit=100)
    return r[0]*xs

def nu_s(gamma_e, B):
    nuc = gamma_e*gamma_e*nu_L(B)
    return nuc

def nu_c(gamma_e, B):
    nuc = (3./2)*gamma_e*gamma_e*nu_L(B)
    return nuc

def get_nu_i_min(gamma_e, B, gamma_min, gamma, fre):
    nus = nu_s(gamma_e, B)
    nus_min = 1.2*1e6*np.power(gamma_min, 2)*B
    f_nu_i_min = []
    for elm in fre:
        f_nu_i_min.append(
            nus_min*(nus)/(4*np.power(gamma, 2)*(1-(h*elm)/(gamma*mec2))))
    max_f_nu_i_min = max(f_nu_i_min)
    return max_f_nu_i_min

def get_nu_i_max(gamma_e, B, gamma_max, gamma, fre):
    nus = nu_s(gamma_e, B)
    nus_max = 1.2*1e6*np.power(gamma_max, 2)*B
    f_nu_i_max = []
    for elm in fre:
        f_nu_i_max.append(nus_max*(nus)/(1-(h*elm)/(gamma*mec2)))
    max_f_nu_i_max = max(f_nu_i_max)
    return max_f_nu_i_max

def luminosity_ssc(R, j_ssc):
    V_R = (4/3)*np.pi*np.power(R, 3)
    L_ssc = (4*np.pi*V_R*j_ssc)
    return L_ssc

def singleElecSynchroPower(nu, gamma, B):
    nuL = nu_L(B)
    n1 = 2.*np.pi*np.sqrt(3.)*e*e*nuL/c
    nus = nu_c(gamma, B)
    x0 = nu/nus
    y0 = bessel(x0)
    P = y0
    return P

def syncEmissivityKernPL(gamma, nu, p, B):
    ne = np.power(gamma, -p)
    k1 = ne*singleElecSynchroPower(nu, gamma, B)
    return k1

def F_syn(nu_list_SYN, B, p, R, n0, gamma_min, gamma_max, d, dl):

    gamma_list = getLogGammaArray(np.log10(gamma_min), np.log10(gamma_max), 100)
    nuL = nu_L(B)
    coeff_P_syn = (2.*np.pi*np.sqrt(3.)*e*e*nuL)/c
    k0 = (p+2)/(8*np.pi*me)
    flux_syn = []
    freq_plot = []
    assorb = []
    V_R = (4/3)*np.pi*np.power(R, 3)
    for elm_nu in nu_list_SYN: 
        func1 = []
        func2 = []
        for elm_gamma in gamma_list: 
            # x_1 = elm_nu/nu_c(elm_gamma, B)
            # alpha_PL = ((3*sigmaT)/(64*np.pi*me))*(np.power(nuL, (p-2)/2))*(np.power(elm_nu, -(p+4)/2))
            # tau = alpha_PL*R
            asso = singleElecSynchroPower(elm_nu, elm_gamma, B)*pow(elm_gamma, -(p+1))
            js = syncEmissivityKernPL(elm_gamma, elm_nu, p, B)
            func1.append(js)
            func2.append(asso)
        integral_simpson = integrate.simps(func1)
        al = integrate.simps(func2)
        if integral_simpson > 1e-50:
            I = coeff_P_syn*integral_simpson*n0 # (r[0]/alpha)*(1.-exp(-tau))
            # alpha_PL = ((3*sigmaT)/(64*np.pi*me))*(np.power(nuL, (p-2)/2))*(np.power(elm_nu, -(p+4)/2))
            alpha = coeff_P_syn*al*k0*n0/pow(elm_nu, 2.)
            tau = alpha*R
            if tau > 0.1:
                assorb.append(elm_nu*I*R)
                I = I*R*(1.-np.exp(-tau))/tau
                flux_syn.append(elm_nu*I)
                freq_plot.append(elm_nu)
            else:
                flux_syn.append(elm_nu*I*R)
                # flux_syn.append((elm_nu*I*np.power(d,4))/(4*np.pi*np.power(dl,2)))
                assorb.append(elm_nu*I*R)
                freq_plot.append(elm_nu)
    # freq_plot_syn = np.log10(freq_plot)
    # flux_plot_syn = np.log10(flux_syn)
    # plt.plot(np.log10(freq_plot), np.log10(flux_syn), c='blue')
    # plt.ylim(-10, 10)
    # plt.xlim(7., 25.)
    # plt.show()
    return np.log10(freq_plot), np.log10(flux_syn)

def GAMMA_fc(gamma):
    f1 = lambda nui: (4*gamma*h*nui)/(mec2)
    return f1

def q_fc(nu, gamma):
    f0 = lambda nui: (nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))
    return f0

def A1(nu, gamma): 
    a1 = 2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2))
    return a1

def P_ssc(B, nu, p, gamma_min, gamma_max, gamma):

    # asso, tau, js = F_syn(nu_list_SYN, B, p, R, n0, gamma_min, gamma_max)
    # syn_part = (1-(3/4)*tau)*asso

    A = 8*np.pi*np.power(re,2)*c*h
    nu_s_min = 1.2*1e6*np.power(gamma_min, 2)*B
    nu_s_max = 1.2*1e6*np.power(gamma_max, 2)*B
    nus = nu_s(gamma, B)
    nui_min = max(nu_s_min, (nus)/(4*np.power(gamma, 2)*(1-(h*nu)/(gamma*mec2))))
    nui_max = max(nu_s_max, (nus)/(1-(h*nu)/(gamma*mec2)))
    
    '''
    f2 = lambda nui: syncEmissivityKernPL(gamma, nui, p, B)
    f0 = lambda nui: (2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2)))*2*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*np.log((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*f2(nui)*(1/np.power(nui,2))
    f1 = lambda nui: (2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2)))*f2(nui)*(1/np.power(nui,2))
    f3 = lambda nui: (2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2)))*f2(nui)*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*(1/np.power(nui,2))
    f4 = lambda nui: (2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2)))*-2*f2(nui)*(np.power((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))),2))*(1/np.power(nui,2))
    f5 = lambda nui: (np.pi*np.power(re,2)*c*(nu/np.power(gamma,2)))*f2(nui)*(np.power((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))),2))*(np.power((4*gamma*h*nui)/(mec2),2))*(1-(nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*(1/(1-((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*((4*gamma*h*nui)/(mec2))))*(1/np.power(nui,2))
    '''

    f2 = lambda nui: syncEmissivityKernPL(gamma, nui, p, B)
    f0 = lambda nui: f2(nui)*2*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*np.log((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*(1/np.power(nui,2))
    f1 = lambda nui: f2(nui)*(1/np.power(nui,2))
    f3 = lambda nui: f2(nui)*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*(1/np.power(nui,2))
    f4 = lambda nui: f2(nui)*-2*(np.power((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))),2))*(1/np.power(nui,2))
    f5 = lambda nui: f2(nui)*(np.power((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))),2))*(np.power((4*gamma*h*nui)/(mec2),2))*\
        (1-(nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))*(1/(1-((nu)/(4*np.power(gamma, 2)*\
        nui*(1-(h*nu)/(gamma*mec2))))*((4*gamma*h*nui)/(mec2))))*(1/np.power(nui,2))

    f_tot = lambda nui: (f0(nui) + f1(nui) + f3(nui) + f4(nui) + f5(nui))

    r = integrate.quad(f_tot, nui_min, nui_max)

    return r[0]

def emissivity_SSC(B, nu, p, gamma_min, gamma_max, gamma):
    ne = np.power(gamma, -p)
    #q_nui = lambda nui: (nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))
    r1 = (2*np.pi*np.power(re,2)*c)*ne*P_ssc(B, nu, p, gamma_min, gamma_max, gamma)
    # (2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2)))
    # print('j_ssc:', r1)
    return r1

'''
def emissivity_SSC_ghisellini(p, R, n0, gamma, nu): 
    alpha = (p-1)/2
    tau = sigmaT*R*n0
    emissivity_syn = (np.power((4/3), alpha))*(np.power((4/3), alpha-1)/2)*np.power(tau, 2)*(1/(4*np.pi))*(1/4)*(c/R)*np.power(4/3, -alpha)*np.power(gamma, -alpha)*np.power(nu, -alpha)*np.log(1e7/1e15)*np.log(1e13, 1e20)
    return emissivity_syn
'''

def sympyintegral(B, R, d, dl, gamma_min, gamma_max, n0, p, nu_list_SSC):
    
    flux_ssc = []
    freq_plot = []
    gamma_list = np.linspace(gamma_min, gamma_max, 200)
    V_R = (4/3)*np.pi*np.power(R, 3)

    for elm_nu in tqdm(nu_list_SSC):
        j_ssc_list = [] 
        for elm_gamma in tqdm(gamma_list):
            j_ssc = emissivity_SSC(B, elm_nu, p, gamma_min, gamma_max, elm_gamma)
            # print('j_ssc', j_ssc)
            j_ssc_list.append(j_ssc)
            # for j in j_ssc_list: 
                # if np.isnan(j) == True: 
                    # j_ssc_list.remove(j)
        integral_simpson1 = integrate.simps(j_ssc_list)
        # print('integral_simps_ssc_1:', integral_simpson1)
        if np.isnan(integral_simpson1) == True:
            print('integral_simps_ssc_1:', integral_simpson1)
            print('gamma', elm_gamma)
            print('nu:', elm_nu)
            print('j_ssc_list', j_ssc_list)
            # sys.exit()
        if integral_simpson1 > 1e-50:
            # print('integral_simps_ssc_2:', integral_simpson1)
            I1 = integral_simpson1*n0
            print(I1)
            # L_ssc = (4*np.pi*V_R*I1)
            # flux_ssc.append((elm_nu*L_ssc*np.power(d,4))/(4*np.pi*np.power(dl,2)))
            flux_ssc.append(elm_nu*I1*R)
            freq_plot.append(elm_nu)
    # freq_plot_ssc = np.log10(freq_plot)
    # flux_plot_ssc = np.log10(flux_ssc)       
    # plt.plot(np.log10(freq_plot), np.log10(flux_ssc), c='red')
    # plt.xlim(12., 24.)
    # plt.show()
    return np.log10(freq_plot), np.log10(flux_ssc)

def plot_syn_ssc(freq_syn, freq_ssc, flux_syn, flux_ssc):
    
    plt.plot(freq_syn, flux_syn, c='red')
    print(flux_syn)
    plt.plot(freq_ssc, flux_ssc, c='blue')
    print(flux_ssc)
    plt.ylim(-10, 5)
    # plt.xlim(7., 25.)
    plt.show()

def main():

    gamma_e = 100.
    B = 0.1
    gamma_min = 10.
    gamma_max = 100000.
    gammaL= 100.
    dt = 1e3
    z = 1
    thetaObs = 1./gammaL
    n0 = 10e4
    p = 2.5
    # nu = 1e17

    d = doppler(gammaL, thetaObs)
    dl = luminosityDistance(z)
    nu_list_SYN = getLogFreqArray(5., 22., 100)
    # fre = getLogFreqArray(5., 18., 400)
    nu_list_SSC = getLogFreqArray(15., 25., 100)
    # nu_list_SSC = np.logspace(12., 22., 20, endpoint=True)
    R = regionSize(dt, d, z)
    nu_min = get_nu_i_min(gamma_e, B, gamma_min, gammaL, nu_list_SYN)
    nu_max = get_nu_i_max(gamma_e, B, gamma_max, gammaL, nu_list_SYN)

    # integrale = int_jssc(gamma, nu_min, nu_max)

    # freq_plot_syn, flux_plot_syn = F_syn(nu_list_SYN, B, p, R, n0, gamma_min, gamma_max)
    # print(freq_plot_syn, flux_plot_syn)
    # freq_plot_syn, flux_plot_syn = F_syn(nu_list_SYN, B, p, R, n0, gamma_min, gamma_max, d, dl)
    freq_plot_syn, flux_plot_syn = F_syn(nu_list_SYN, B, p, R, n0, gamma_min, gamma_max, d, dl)
    freq_plot_ssc, flux_plot_ssc = sympyintegral(B, R, d, dl, gamma_min, gamma_max, n0, p, nu_list_SSC)
    # F_syn(nu_list_SYN, B, p, R, n0, gamma_min, gamma_max, d, dl)
    # sympyintegral(B, R, d, dl, gamma_min, gamma_max, n0, p, nu_list_SSC)
    plot_syn_ssc(freq_plot_syn, freq_plot_ssc, flux_plot_syn, flux_plot_ssc)

    # F_syn = syn_sympyintegral(fre_syn, Lsyn, d, dl)
    # print(F_syn)
    # P_ssc, j_ssc_ b1, j_ssc_2, j_ssc_3, j_ssc_tot = sympyintegral(A, R, nu_min, nu_max, gamma_min, gamma_max, n0, p)
    # print('j_ssc_tot', j_ssc_tot)

    # L_ssc = luminosity_ssc(R, j_ssc_tot)
    # F_ssc = flux_ssc(fre, L_ssc, d, dl)

    # print(F_ssc)
    # print(type(F_ssc))

    # f_c = Jones_kernel(100)
    # P_ssc = power_ssc(gamma_e, B, gamma_min, gamma_max, gamma)
    # print(P_ssc)
    # j_ssc = emissivity_ssc(n0, gamma_min, gamma_max, p, P_ssc)
    # L_ssc = luminosity_ssc(R, j_ssc)
    # Flux_ssc = flux_ssc(fre, L_ssc, d, dl)


if __name__ == "__main__":
    main()

'''
def I_2(nu_list, B, p, R, n0, gamma_min, gamma_max, nu, gamma, nui):
    tau, js = F_syn(nu_list, B, p, R, n0, gamma_min, gamma_max)
    syn_part = (3./4.)*(R/c)*(1-(3/4)*tau)*js 
    i2 = (nu/4*np.power(gamma,2)*nui)*(1/(h*nui))*((1)+(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))-2*\
    (np.exp(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)**2)))*(syn_part)
    return i2

def I_3(nu_list, B, p, R, n0, gamma_min, gamma_max, nu, gamma, nui): 
    tau, js = F_syn(nu_list, B, p, R, n0, gamma_min, gamma_max)
    syn_part = (3./4.)*(R/c)*(1-(3/4)*tau)*js 
    i3 = (nu/4*np.power(gamma,2)*nui)*(1/(h*nui))*((4*gamma*nui*h/mec2)**2)*(np.exp((nu/4*np.power(gamma,2)*nui\
    *(1-h*nu/gamma*mec2))**2))*(1/(2+2*(4*gamma*nui*h/mec2)*(nu/(4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))))\
    *(1-(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))*(syn_part)
    return i3
'''

# ((f2(nui))/(np.power(nui,2)))*(np.power(((4*gamma*h*nui)/(mec2))*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))),2)*(1-((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))))/(1-((4*gamma*h*nui)/(mec2))*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))))
# (1/(np.power(nui,2)))*(np.power(((4*gamma*h*nui)/(mec2))*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))),2)*(1-((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2))))))/(1-((4*gamma*h*nui)/(mec2))*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))))
'''
I1 = a1*(2*q*np.log(q)*syn_part)/(np.power(nui,2))
I2 = a1*(syn_part)/(np.power(nui,2))
I3 = a1*(syn_part*q)/(np.power(nui,2))
I4 = a1*(-syn_part*2*np.power(q,2))/(np.power(nui,2))
I5 = (a1/2)*((syn_part)/(np.power(nui,2)))*(np.power(Gamma*q,2)*(1-q))/(1-Gamma*q)
'''
# i1 = (nu/4*np.power(gamma,2)*nui)*(2*nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))*(np.log(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))*(1/(h*nui))*(syn_part)
# i2 = (nu/4*np.power(gamma,2)*nui)*(1/(h*nui))*((1)+(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))-2*(np.exp(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)**2)))*(syn_part)
# i3 = (nu/4*np.power(gamma,2)*nui)*(1/(h*nui))*((4*gamma*nui*h/mec2)**2)*(np.exp((nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2))**2))*(1/(2+2*(4*gamma*nui*h/mec2)*(nu/(4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))))*(1-(nu/4*np.power(gamma,2)*nui*(1-h*nu/gamma*mec2)))*(syn_part)
'''
f_f2 = lambda nui: 2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2))*(f2(nui)*(2*(nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))*np.log((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))))/(np.power(nui,2)))+\
    (f2(nui))/(np.power(nui,2))+f2(nui)*(1*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))))/(np.power(nui,2))+\
    f2(nui)*(-1*2*np.power(((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))),2))/(np.power(nui,2))
f = lambda nui: 2*np.pi*np.power(re,2)*c*(nu/np.power(gamma,2))*(1*(2*(nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))*np.log((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))))/(np.power(nui,2)))+\
    (1)/(np.power(nui,2))+1*(1*((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))))/(np.power(nui,2))+\
    (-1*2*np.power(((nu)/(4*np.power(gamma, 2)*nui*(1-(h*nu)/(gamma*mec2)))),2))/(np.power(nui,2))
'''
# a1 = A1(nu, gamma)
# Gamma = GAMMA_fc(gamma)
# print(Gamma)
# q = q_fc(nu, gamma)
# print(q)
# f = lambda nui: a1*(2*q*np.log(q)*syn_part)/(np.power(nui,2))+a1*(syn_part)/(np.power(nui,2))+a1*(syn_part*q)/(np.power(nui,2))+\
    # a1*(-syn_part*2*np.power(q,2))/(np.power(nui,2))+(a1/2)*((syn_part)/(np.power(nui,2)))*(np.power(Gamma*q,2)*(1-q))/(1-Gamma*q)
# r = integrate.quad(f, nu_min, nu_max)
# f = lambda nui: a1*(2*q*np.log(q))/(np.power(nui,2))+a1*(1)/(np.power(nui,2))+a1*(1*q)/(np.power(nui,2))+\
    # a1*(-1*2*np.power(q,2))/(np.power(nui,2))+(a1/2)*((1)/(np.power(nui,2)))*(np.power(Gamma*q,2)*(1-q))/(1-Gamma*q)

