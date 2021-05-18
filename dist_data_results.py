import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
import sys


path_sim = '/Users/alessandraberretta/JetFit/data_100_theta_obs_0.017_20/summary_P1/'
path_fit = '/Users/alessandraberretta/JetFit/results_100_theta_obs_P1_0.017_20/'


E_f = []
E_s = []
err_pos_E_f = []
err_neg_E_f = []
Eta0_f = []
Eta0_s = []
err_pos_Eta0_f = []
err_neg_Eta0_f = []
GammaB_f = []
GammaB_s = []
err_pos_GammaB_f = []
err_neg_GammaB_f = []
epsB_f = []
epsB_s = []
err_pos_epsB_f = []
err_neg_epsB_f = []
epse_f = []
epse_s = []
err_pos_epse_f = []
err_neg_epse_f = []
n_f = []
n_s = []
err_pos_n_f = []
err_neg_n_f = []
p_f = []
p_s = []
err_pos_p_f = []
err_neg_p_f = []
theta_obs_f = []
theta_obs_s = []
err_pos_theta_obs_f = []
err_neg_theta_obs_f = []
p_s_low = []
# E_dev = []


list_file_sim = [path_sim +
                 file for file in os.listdir(path_sim) if file.startswith('summary')]
list_file_fit = [path_fit +
                 file for file in os.listdir(path_fit) if file.startswith('summary')]


for elm in list_file_fit:

    # num_P = elm[elm.rfind('P')+1:elm.rfind('_')]
    num_P = elm[elm.rfind('/')+1:].split('_')[2]
    # file_sim = [s for s in list_file_sim if f"{num_P}" in s and s[s.rfind('P'):s.rfind('_E')] == num_P]
    file_sim = [s for s in list_file_sim if f"{num_P}" in s and s[s.rfind(
        '/')+1:].split('_')[1] == num_P]

    df = pd.read_csv(file_sim[0], sep='\t')
    E_s.append(df['Values'][3])
    Eta0_s.append(df['Values'][0])
    GammaB_s.append(df['Values'][1])
    epsB_s.append(df['Values'][4])
    epse_s.append(df['Values'][5])
    n_s.append(df['Values'][6])
    p_s.append(df['Values'][7])
    theta_obs_s.append(df['Values'][2])

    df = pd.read_csv(elm, sep='\t')
    '''
    E_f.append(df['Values'][0])
    err_pos_E_f.append(df['err_pos'][0])
    err_neg_E_f.append(df['err_neg'][0])
    Eta0_f.append(df['Values'][0])
    err_pos_Eta0_f.append(df['err_pos'][0])
    err_neg_Eta0_f.append(df['err_neg'][0])
    GammaB_f.append(df['Values'][1])
    err_pos_GammaB_f.append(df['err_pos'][1])
    err_neg_GammaB_f.append(df['err_neg'][1])
    epsB_f.append(df['Values'][0])
    err_pos_epsB_f.append(df['err_pos'][0])
    err_neg_epsB_f.append(df['err_neg'][0])
    epse_f.append(df['Values'][0])
    err_pos_epse_f.append(df['err_pos'][0])
    err_neg_epse_f.append(df['err_neg'][0])
    n_f.append(df['Values'][1])
    err_pos_n_f.append(df['err_pos'][1])
    err_neg_n_f.append(df['err_neg'][1])
    p_f.append(df['Values'][0])
    err_pos_p_f.append(df['err_pos'][0])
    err_neg_p_f.append(df['err_neg'][0])
    '''
    theta_obs_f.append(df['Values'][0])
    err_pos_theta_obs_f.append(df['err_pos'][0])
    err_neg_theta_obs_f.append(df['err_neg'][0])
'''
    
for i in range(len(E_s)):
    if np.log10(E_s[i]) > 0.06 and np.log10(E_s[i]) < 0.065 and E_f[i] > 0.06:
        # if np.log10(E_s[i]) - E_f[i] < 0.001:
        print(list_file_fit[i])

for i in range(len(epsB_s)):
    if np.log10(epsB_s[i]) < -3.03 and np.log10(epsB_s[i]) > -3.10 and epsB_f[i] < -3.6 and epsB_f[i] > -3.8:
        # if np.log10(E_s[i]) - E_f[i] < 0.001:
        print(list_file_fit[i])


p_low = []
for idx, val in enumerate(p_s):
    if val < 2.15:
        p_low.append(val)

for elm in list_file_sim:
    df = pd.read_csv(elm, sep='\t')
    if df['Values'][7] < 2.15:
        p_s_low.append(elm)

# print(p_low, len(p_low), p_s_low, len(p_s_low))


E_f_pow = np.power(10, E_f)
E_s_pow = np.power(10, E_s)
err_pos_E_f_pow = np.power(10, err_pos_E_f)
err_neg_E_f_pow = np.power(10, err_neg_E_f)
epsB_f_pow = np.power(10, epsB_f)
epsB_s_pow = np.power(10, epsB_s)
err_pos_epsB_f_pow = np.power(10, err_pos_epsB_f)
err_neg_epsB_f_pow = np.power(10, err_neg_epsB_f)
epse_f_pow = np.power(10, epse_f)
epse_s_pow = np.power(10, epse_s)
err_pos_epse_f_pow = np.power(10, err_pos_epse_f)
err_neg_epse_f_pow = np.power(10, err_neg_epse_f)
n_f_pow = np.power(10, n_f)
n_s_pow = np.power(10, n_s)
err_pos_n_f_pow = np.power(10, err_pos_n_f)
err_neg_n_f_pow = np.power(10, err_neg_n_f)


err_E_f = (np.array(err_pos_E_f) - np.array(err_neg_E_f))/2
err_E_s = 0.13379407681190938*np.ones(len(err_E_f))
err_p_f = (np.array(err_pos_p_f) - np.array(err_neg_p_f))/2
err_p_s = 0.4195430461130414*np.ones(len(err_p_f))
# err_p_s = 0.2*np.array(p_s)
err_Eta0_f = abs((np.array(err_pos_Eta0_f) + np.array(err_neg_Eta0_f))/2)
err_Eta0_s = 0.4603489203236277*np.ones(len(err_Eta0_f))
err_GammaB_f = abs((np.array(err_pos_GammaB_f) + np.array(err_neg_GammaB_f))/2)
err_GammaB_s = 0.57*np.ones(len(err_GammaB_f))
err_epsb_f = abs((np.array(err_pos_epsB_f) + np.array(err_neg_epsB_f))/2)
err_epsb_s = 0.00015384891068757133*np.ones(len(err_epsb_f))
err_epse_f = (np.array(err_pos_epse_f) - np.array(err_neg_epse_f))/2
err_epse_s = 0.004702845692750316*np.ones(len(err_epse_f))
err_n_f = abs((np.array(err_pos_n_f) + np.array(err_neg_n_f))/2)
err_n_s = 0.05699152118*np.asarray(n_s)
err_theta_obs_f = abs((np.array(err_pos_theta_obs_f) +
                       np.array(err_neg_theta_obs_f))/2)
err_theta_obs_s = 0.00953959783*np.ones(len(err_theta_obs_f))
'''

fig = plt.figure()
ax = plt.gca()
# fig, ax = plt.subplots()

'''
gamma_f = []
for num1, num2 in zip(Eta0_f, GammaB_f):
    gamma_f.append(2 * num1 * num2)

gamma_s = []
for num3, num4 in zip(Eta0_s, GammaB_s):
    gamma_s.append(2 * num3 * num4)


gamma = np.asarray(Eta0_f)*np.asarray(GammaB_f)


# scatterplots
E_dev = (np.log10(E_s) - np.asarray(E_f))/np.asarray(E_s)
plt.hist(E_dev, bins=20, histtype='step', color='black')
ax.set_xlabel(
    r'$\frac{E\_s - E\_f}{E\_s}$', fontsize='x-large', alpha=1)
mean = np.mean(E_dev)
sigma = np.std(E_dev)
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean, ),
    r'$\sigma=%.2f$' % (sigma, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='black')
plt.show()
# plt.scatter(E_s, np.power(10, E_f), c='black')
ax.errorbar(np.log10(E_s), E_f, yerr=err_pos_E_f,
            xerr=None, fmt='o', ms=3.5, mec='royalblue', c='black')
plt.xlabel('E_s')
plt.ylabel('E_f')
# ax.set_xlim(-0.50, 0.50)
# ax.set_ylim(-0.50, 0.50)
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='royalblue')
print(ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='royalblue'))
plt.title('E comparison')
# plt.savefig('/Users/alessandraberretta/JetFit/plot/E_eb_P1_0.15_40.pdf')

plt.scatter(Eta0_s, Eta0_f, c='black')
plt.xlabel('Eta0_s')
plt.ylabel('Eta0_f')
plt.title('Eta0 comparison')
plt.savefig('/Users/alessandraberretta/JetFit/plot/Eta0_scatter_P1_1par.pdf')

plt.scatter(GammaB_s, GammaB_f, c='black')
plt.xlabel('GammaB_s')
plt.ylabel('GammaB_f')
plt.title('GammaB comparison')
plt.savefig('/Users/alessandraberretta/JetFit/plot/GammaB_scatter_P1_1par.pdf')

epsb_dev = (np.asarray(epsB_s) - np.power(10, epsB_f))/np.asarray(epsB_s)
plt.hist(epsb_dev, bins=30, histtype='step', color='black')
ax.set_xlabel(
    r'$\frac{epsB\_s - epsB\_f}{epsB\_s}$', fontsize='x-large', alpha=1)
mean = np.mean(epsb_dev)
sigma = np.std(epsb_dev)
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean, ),
    r'$\sigma=%.2f$' % (sigma, )))
ax.text(0.04, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='black')
plt.show()
# plt.scatter(epsB_s, np.power(10, epsB_f), c='black')

ax.errorbar(np.log10(epsB_s), epsB_f, yerr=err_pos_epsB_f,
            xerr=None, fmt='o', ms=4.5, mec='royalblue', c='black')
plt.xlabel('epsB_s')
plt.ylabel('epsB_f')
ax.set_xlim(-3.20, -2.60)
ax.set_ylim(-3.20, -2.60)
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='royalblue')
plt.title('epsB comparison')
# plt.savefig('/Users/alessandraberretta/JetFit/plot/epsb_scatter_P1_1par_eb.pdf')
epse_dev = (np.asarray(epse_s) - np.power(10, epse_f))/np.asarray(epse_s)
plt.hist(epse_dev, bins=30, histtype='step', color='black')
ax.set_xlabel(
    r'$\frac{epse\_s - epse\_f}{epse\_s}$', fontsize='x-large', alpha=1)
mean = np.mean(epse_dev)
sigma = np.std(epse_dev)
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean, ),
    r'$\sigma=%.2f$' % (sigma, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='black')
plt.show()
# plt.scatter(epse_s, np.power(10, epse_f), c='black')
ax.errorbar(np.log10(epse_s), epse_f, yerr=err_pos_epse_f,
            xerr=None, fmt='o', ms=4.5, mec='royalblue', c='black')
plt.xlabel('epse_s')
plt.ylabel('epse_f')
ax.set_xlim(-0.65, -0.25)
ax.set_ylim(-0.65, -0.25)
ax.plot([0, 1], [0, 1], transform=ax.transAxes)
plt.title('epse comparison')
plt.savefig('/Users/alessandraberretta/JetFit/plot/epse_eb_P1_0.4_20.pdf')

plt.scatter(n_s, np.power(10, n_f), c='black')
# ax.set_xlim(0.11, 1.20)
# ax.set_ylim(0.11, 1.20)
plt.xlabel('n_s')
plt.ylabel('n_f')
plt.title('n comparison')
# plt.savefig('/Users/alessandraberretta/JetFit/plot/n_scatter_P1_1par.pdf')

# plt.scatter(p_s, p_f, c='black')
ax.errorbar(p_s, p_f, yerr=err_pos_p_f,
            xerr=None, fmt='o', ms=3.5, mec='royalblue', c='black')
plt.xlabel('p_s')
plt.ylabel('p_f')
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='royalblue')
# ax.set_xlim(2.8, 4.6)
# ax.set_ylim(2.8, 4.6)
# ax.set_xlim(2.1, 4)
# ax.set_ylim(2.1, 4)
plt.title('p comparison')
# plt.savefig('/Users/alessandraberretta/JetFit/plot/p_scatter_P1_all_3_4_diag.pdf')
'''
plt.scatter(theta_obs_s, theta_obs_f, c='black')
plt.xlabel('theta_obs_s')
plt.ylabel('theta_obs_f')
# ax.set_xlim(0.06, 0.23)
# ax.set_ylim(0.06, 0.23)
plt.title('theta_obs comparison')
# plt.savefig('/Users/alessandraberretta/JetFit/plot/theta_obs_scatter_P1_1par.pdf')
'''
plt.scatter(gamma_s, gamma_f, c='black')
plt.xlabel('gamma_s')
plt.ylabel('gamma_f')
ax.set_xlim(10, 200)
ax.set_ylim(10, 200)
plt.title('gamma_tot comparison')
# plt.savefig('/Users/alessandraberretta/JetFit/plot/gamma_scatter_P1_1par_2.pdf')


# hisotgrams

# print(np.min(E_f), np.max(E_f))
# plt.hist(np.array(E_f), bins=50, histtype='step', color='red')
# plt.hist(np.log10(E_s), bins=50, histtype='step', color='blue')
plt.hist(err_pos_E_f, bins=50, histtype='step', color='black')
ax.set_xlabel('Energy (/10^50 erg)')
mean = np.mean(err_pos_E_f)
sigma = np.std(err_pos_E_f)
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean, ),
    r'$\sigma=%.2f$' % (sigma, )))
ax.text(0.65, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='black')
# ax.text(0.65, 0.95, 'fitted', transform=ax.transAxes, fontsize=15, verticalalignment='top', color='red')
# ax.text(0.65, 0.90, 'simulated', transform=ax.transAxes,fontsize=15, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/E_hist.pdf')

plt.hist(Eta0_f, bins=60, histtype='step', color='red')
plt.hist(Eta0_s, bins=60, histtype='step', color='blue')
ax.set_xlabel('Eta0')
ax.text(0.05, 0.95, 'fitted', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.05, 0.90, 'simulated', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/Eta0_hist.pdf')

plt.hist(GammaB_f, bins=80, histtype='step', color='red')
plt.hist(GammaB_s, bins=80, histtype='step', color='blue')
ax.set_xlabel('GammaB')
ax.text(0.05, 0.95, 'fitted', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.05, 0.90, 'simulated', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/GammB_hist.pdf')

plt.hist(epsB_f, bins=40, histtype='step', color='red')
plt.hist(np.log10(epsB_s), bins=30, histtype='step', color='blue')
ax.set_xlabel('epsB')
# ax.set_xscale('log')
ax.text(0.05, 0.95, 'fitted', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.05, 0.90, 'simulated', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
plt.savefig(
    '/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/epsb_hist.pdf')
plt.show()

plt.hist(epse_f, bins=20, histtype='step', color='red')
plt.hist(np.log10(epse_s), bins=20, histtype='step', color='blue')
ax.set_xlabel('epse')
# ax.set_xscale('log')
ax.text(0.05, 0.95, 'fitted', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.05, 0.90, 'simulated', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/epse_hist.pdf')

plt.hist(n_f, bins=30, histtype='step', color='red')
plt.hist(np.log10(n_s), bins=30, histtype='step', color='blue')
ax.set_xlabel('n')
ax.text(0.05, 0.95, 'fitted', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.05, 0.90, 'simulated', transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
plt.savefig(
    '/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/n_hist.pdf')

plt.hist(p_low, bins=10, histtype='step', color='red')
# plt.hist(p_s, bins=20, histtype='step', color='blue')
ax.set_xlabel('p')
# plt.savefig('p.png')

plt.hist(theta_obs_f, bins=30, histtype='step', color='red')
plt.hist(theta_obs_s, bins=30, histtype='step', color='blue')
ax.set_xlabel('theta_obs')


# discrepancy plots
a = np.power(10, E_f) - np.array(E_s)
b = np.sqrt(pow(np.power(10, err_E_f), 2) + pow(err_E_s, 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 40
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlabel(
    r'$\frac{E\_f - E\_s}{\sqrt{(err\_E\_f)^{2} + (err\_E\_s)^{2}}}$', fontsize='x-large', alpha=1)
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.65, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.65, 0.80, textstr_t, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/E_discrepancy.pdf')

a = np.asarray(Eta0_f) - np.asarray(Eta0_s)
b = np.sqrt(pow(np.asarray(err_Eta0_f), 2) + pow(np.asarray(err_Eta0_s), 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
err_mean_s = sigma_s/math.sqrt(len(r))
num_bins = 40
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
extension_l = np.linspace(-6, -4.2, 20)
extension_r = np.linspace(2.2, 6, 20)
bins_tot = list(extension_l) + list(bins) + list(extension_r)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (np.asarray(bins_tot) - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (np.asarray(bins_tot) - mean_s))**2))
ax.plot(np.asarray(bins_tot), y_t, color='blue')
ax.plot(np.asarray(bins_tot), y_s, color='red')
ax.set_xlim(-8, 8)
ax.set_xlabel(
    r'$\frac{Eta0\_f - Eta0\_s}{\sqrt{(err\_Eta0\_f)^{2} + (err\_Eta0\_s)^{2}}}$', fontsize='x-large', alpha=1)
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, ), r'$err\_mean\_s=%.2f$' % (err_mean_s, )))
ax.text(0.01, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=12, verticalalignment='top', color='red')
ax.text(0.01, 0.75, textstr_t, transform=ax.transAxes,
        fontsize=12, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/Eta0_discrepancy.png')

a = np.asarray(GammaB_f) - np.asarray(GammaB_s)
b = np.sqrt(pow(np.asarray(err_GammaB_f), 2) +
            pow(np.asarray(err_GammaB_s), 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
err_mean_s = sigma_s/math.sqrt(len(r))
num_bins = 30
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlim(-6, 6)
ax.set_xlabel(
    r'$\frac{GammaB\_f - GammaB\_s}{\sqrt{(err\_GammaB\_f)^{2} + (err\_GammaB\_s)^{2}}}$', fontsize='x-large', alpha=2)
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, ), r'$err\_mean\_s=%.2f$' % (err_mean_s, )))
ax.text(0.02, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.02, 0.78, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/GammaB_discrepancy.pdf')

a = np.array(p_f) - np.array(p_s)
b = np.sqrt(pow(np.asarray(err_pos_p_f), 2) + pow(np.asarray(err_p_s), 2))
r = [a[i] / b[i] for i in range(len(a))]
r_array = np.asarray(r)
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 40
n, bins, patches = plt.hist(
    r_array, num_bins, histtype='step', color='black', density=True)
# extension_l = np.linspace(-4, -1.29, 20)
# extension_r = np.linspace(2.3, 4, 20)
# bins_tot = list(extension_l) + list(bins) + list(extension_r)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
# ax.set_xlim(-0.2, 0.2)
ax.set_xlabel(
    r'$\frac{p\_f - p\_s}{\sqrt{(err\_p\_f)^{2} + (err\_p\_s)^{2}}}$', fontsize='x-large', alpha=1)
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.05, 0.85, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
# plt.savefig('p_discrepancy_P2.pdf')
plt.scatter(b, a, c='black')

a = np.asarray(n_f) - np.asarray(np.log10(n_s))
b = np.sqrt(pow(np.asarray(err_n_f), 2) +
            pow(np.asarray(np.log10(err_n_s)), 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 40
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlim(-5, 5)
ax.set_xlabel(
    '(n_f - n_s)/sqrt((err_n_f)^2 + (err_n_s)^2)')
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.05, 0.85, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
# plt.savefig('n_discrepancy.pdf')

a = np.asarray(theta_obs_f) - np.asarray(theta_obs_s)
b = np.sqrt(pow(np.asarray(err_theta_obs_f), 2) +
            pow(np.asarray(err_theta_obs_s), 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 30
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlabel(
    r'$\frac{theta\_obs\_f - theta\_obs\_s}{\sqrt{(err\_theta\_obs\_f)^{2} + (err\_theta\_obs\_s)^{2}}}$', fontsize='x-large', alpha=1)
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.05, 0.85, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
# plt.savefig('theta_obs_discrepancy.pdf')

a = np.asarray(epsB_f) - np.asarray(np.log10(epsB_s))
b = np.sqrt(pow(np.asarray(err_epsb_f), 2) +
            pow(np.asarray(np.log10(err_epsb_s)), 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(a)
sigma_s = np.std(a)
num_bins = 40
n, bins, patches = plt.hist(
    a, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlabel(
    '(epsB_f - epsB_s)/|err_pos_epsB_f - err_neg_epsB_f|')
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.05, 0.85, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/epsb_discrepancy.pdf')

a = np.asarray(epse_f) - np.asarray(np.log10(epse_s))
b = np.sqrt(pow(np.array(err_epse_f), 2) +
            pow(np.array(np.log10(err_epse_s)), 2))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 30
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
# ax.set_xlabel('(epse_f - epse_s)/|err_pos_epse_f - err_neg_epse_f|')
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.65, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.65, 0.85, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/epse_discrepancy.pdf')


# histogram + scatter plots
def scatter_hist(x, y, ax, ax_histx, ax_histy):

    ax_histx.tick_params(axis="x", labelbottom=True)
    ax_histy.tick_params(axis="y", labelleft=True)

    # the scatter plot:
    ax.scatter(x, y, color='black')
    # ax.set_xlim(-1, 1)
    # ax.set_ylim(-2, 2)
    ax.set_xlabel('epsb_s')
    ax.set_ylabel('epsb_f')

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    bins_x = 40
    bins_y = 40

    ax_histx.hist(x, bins=bins_x, histtype='step', color='black')
    ax_histy.hist(y, bins=bins_y, orientation='horizontal',
                  histtype='step', color='black')
    # ax.set_xscale('log')


# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes(rect_scatter)

ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# scatter_hist(np.log10(E_s), E_f, ax, ax_histx, ax_histy)
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/E_scatter_dist.pdf')

# scatter_hist(Eta0_s, Eta0_f, ax, ax_histx, ax_histy)
# plt.savefig('/Users/alessandraberretta/JetFit/Eta0_scatter_dist_2.pdf')

# scatter_hist(GammaB_s, GammaB_f, ax, ax_histx, ax_histy)
# plt.savefig('/Users/alessandraberretta/JetFit/GammaB_scatter_dist.pdf')

# scatter_hist(p_s, p_f, ax, ax_histx, ax_histy)
# plt.savefig('p_scatter_dist_P2.pdf')

# scatter_hist(np.log10(n_s), n_f, ax, ax_histx, ax_histy)
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/n_scatter_dist.pdf')

# scatter_hist(theta_obs_s, theta_obs_f, ax, ax_histx, ax_histy)

scatter_hist(np.log10(epsB_s), epsB_f, ax, ax_histx, ax_histy)
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/epsb_scatter_dist.pdf')

# scatter_hist(epse_s, epse_f, ax, ax_histx, ax_histy)
# plt.savefig('/Users/alessandraberretta/JetFit/plot_test_analysis_250lc_200/epse_scatter_dist.pdf')
'''
plt.show()
