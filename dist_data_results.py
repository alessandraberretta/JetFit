import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
import sys


path_sim = '/Users/alessandraberretta/JetFit/data_200_lc/summary/'
path_fit = '/Users/alessandraberretta/JetFit/data_results_100/'


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


list_file_sim = [path_sim +
                 file for file in os.listdir(path_sim) if file.startswith('summary')]
list_file_fit = [path_fit +
                 file for file in os.listdir(path_fit) if file.startswith('summary')]

for elm in list_file_fit:

    num_P = elm[elm.rfind('P')+1:elm.rfind('_')]

    file_sim = [s for s in list_file_sim if f"P{num_P}" in s and s[s.rfind(
        'P')+1:s.rfind('.')] == num_P]

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
    E_f.append(df['Values'][3])
    err_pos_E_f.append(df['err_pos'][3])
    err_neg_E_f.append(df['err_neg'][3])
    Eta0_f.append(df['Values'][0])
    err_pos_Eta0_f.append(df['err_pos'][0])
    err_neg_Eta0_f.append(df['err_neg'][0])
    GammaB_f.append(df['Values'][1])
    err_pos_GammaB_f.append(df['err_pos'][1])
    err_neg_GammaB_f.append(df['err_neg'][1])
    epsB_f.append(df['Values'][4])
    err_pos_epsB_f.append(df['err_pos'][4])
    err_neg_epsB_f.append(df['err_neg'][4])
    epse_f.append(df['Values'][5])
    err_pos_epse_f.append(df['err_pos'][5])
    err_neg_epse_f.append(df['err_neg'][5])
    n_f.append(df['Values'][6])
    err_pos_n_f.append(df['err_pos'][6])
    err_neg_n_f.append(df['err_neg'][6])
    p_f.append(df['Values'][7])
    err_pos_p_f.append(df['err_pos'][7])
    err_neg_p_f.append(df['err_neg'][7])
    theta_obs_f.append(df['Values'][2])
    err_pos_theta_obs_f.append(df['err_pos'][2])
    err_neg_theta_obs_f.append(df['err_neg'][2])


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


fig = plt.figure()
ax = plt.gca()


# scatterplots
'''
plt.scatter(E_s_pow, E_f_pow, c='black')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('E_s')
plt.ylabel('E_f')
plt.title('E comparison')

plt.scatter(Eta0_s, Eta0_f, c='black')
# ax.set_yscale('log')
# ax.set_xscale('log')
plt.xlabel('Eta0_s')
plt.ylabel('Eta0_f')
plt.title('Eta0 comparison')

plt.scatter(GammaB_s, GammaB_f, c='black')
# ax.set_yscale('log')
# ax.set_xscale('log')
plt.xlabel('GammaB_s')
plt.ylabel('GammaB_f')
plt.title('GammaB comparison')

plt.scatter(epsB_s_pow, epsB_f_pow, c='black')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('epsB_s')
plt.ylabel('epsB_f')
plt.title('epsB comparison')

plt.scatter(epse_s_pow, epse_f_pow, c='black')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('epse_s')
plt.ylabel('epse_f')
plt.title('epse comparison')

plt.scatter(n_s_pow, n_f_pow, c='black')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('n_s')
plt.ylabel('n_f')
plt.title('n comparison')

plt.scatter(p_s, p_f, c='black')
# ax.set_yscale('log')
# ax.set_xscale('log')
plt.xlabel('p_s')
plt.ylabel('p_f')
plt.title('p comparison')

plt.scatter(theta_obs_s, theta_obs_f, c='black')
# ax.set_yscale('log')
# ax.set_xscale('log')
plt.xlabel('theta_obs_s')
plt.ylabel('theta_obs_f')
plt.title('theta_obs comparison')



# hisotgrams
plt.hist(E_f_pow, bins=55, histtype='step', color='red')
plt.hist(E_s_pow, bins=15, histtype='step', color='blue')
ax.set_xscale('log')
ax.set_xlabel('Energy (10^50 erg)')
# plt.savefig('E.png')

plt.hist(Eta0_f, bins=15, histtype='step', color='red')
plt.hist(Eta0_s, bins=15, histtype='step', color='blue')
ax.set_xlabel('Eta0')
# plt.savefig('Eta0.png')

plt.hist(GammaB_f, bins=30, histtype='step', color='red')
plt.hist(GammaB_s, bins=30, histtype='step', color='blue')
ax.set_xlabel('GammaB')
# plt.savefig('GammaB.png')

plt.hist(epsB_f_pow, bins=30, histtype='step', color='red')
plt.hist(epsB_s_pow, bins=20, histtype='step', color='blue')
ax.set_xlabel('epsB')
ax.set_xscale('log')
# plt.savefig('epsB.png')

plt.hist(epse_f_pow, bins=25, histtype='step', color='red')
plt.hist(epse_s_pow, bins=30, histtype='step', color='blue')
print(epse_f_pow)
print(epse_s_pow)
ax.set_xlabel('epse')
ax.set_xscale('log')
# plt.savefig('epse.png')

plt.hist(n_f_pow, bins=40, histtype='step', color='red')
plt.hist(n_s, bins=30, histtype='step', color='blue')
print(n_f_pow)
print(n_s)
ax.set_xlabel('n')
ax.set_xscale('log')
# plt.savefig('n.png')

plt.hist(p_f, bins=15, histtype='step', color='red')
plt.hist(p_s, bins=15, histtype='step', color='blue')
ax.set_xlabel('p')
# plt.savefig('p.png')

plt.hist(theta_obs_f, bins=20, histtype='step', color='red')
plt.hist(theta_obs_s, bins=20, histtype='step', color='blue')
ax.set_xlabel('theta_obs')



# discrepancy plots

a = np.asarray(E_f_pow) - np.asarray(E_s_pow)
b = abs(np.asarray(err_pos_E_f_pow) - np.asarray(err_neg_E_f_pow))
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
    '(E_f - E_s)/|err_pos_E_f - err_neg_E_f|')
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.05, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='red')
ax.text(0.05, 0.80, textstr_t, transform=ax.transAxes,
        fontsize=15, verticalalignment='top', color='blue')
plt.xlim(-1.5, 1.5)
plt.savefig('E_discrepancy.pdf')

a = np.asarray(Eta0_f) - np.asarray(Eta0_s)
b = abs(np.asarray(err_pos_Eta0_f) - np.asarray(err_neg_Eta0_f))
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
    '(Eta0_f - Eta0_s)/|err_pos_Eta0_f - err_neg_Eta0_f|')
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
plt.savefig('Eta0_discrepancy.pdf')

a = np.asarray(GammaB_f) - np.asarray(GammaB_s)
b = abs(np.asarray(err_pos_GammaB_f) - np.asarray(err_neg_GammaB_f))
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
    '(GammaB_f - GammaB_s)/|err_pos_GammaB_f - err_neg_GammaB_f|')
textstr_t = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_t, ),
    r'$\sigma=%.2f$' % (sigma_t, )))
textstr_s = '\n'.join((
    r'$\mathrm{mean}=%.2f$' % (mean_s, ),
    r'$\sigma=%.2f$' % (sigma_s, )))
ax.text(0.02, 0.95, textstr_s, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='red')
ax.text(0.02, 0.85, textstr_t, transform=ax.transAxes,
        fontsize=13, verticalalignment='top', color='blue')
plt.savefig('GammaB_discrepancy.pdf')

a = np.asarray(p_f) - np.asarray(p_s)
b = abs(np.asarray(err_pos_p_f) - np.asarray(err_neg_p_f))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 35
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlabel(
    '(p_f - p_s)/|err_pos_p_f - err_neg_p_f|')
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
plt.savefig('p_discrepancy.pdf')

a = np.asarray(n_f_pow) - np.asarray(n_s)
b = abs(np.asarray(err_pos_n_f_pow) - np.asarray(err_neg_n_f_pow))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 50
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlabel(
    '(n_f - n_s)/|err_pos_n_f - err_neg_n_f|')
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
plt.savefig('n_discrepancy.pdf')

a = np.asarray(theta_obs_f) - np.asarray(theta_obs_s)
b = abs(np.asarray(err_pos_theta_obs_f) - np.asarray(err_neg_theta_obs_f))
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
    '(theta_obs_f - theta_obs_s)/|err_pos_theta_obs_f - err_neg_theta_obs_f|')
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
plt.savefig('theta_obs_discrepancy.pdf')

a = np.asarray(epsB_f_pow) - np.asarray(epsB_s_pow)
b = abs(np.asarray(err_pos_epsB_f_pow) - np.asarray(err_neg_epsB_f_pow))
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
plt.savefig('epsB_discrepancy.pdf')
'''
a = np.asarray(epse_f_pow) - np.asarray(epse_s_pow)
b = abs(np.asarray(err_pos_epse_f_pow) - np.asarray(err_neg_epse_f_pow))
r = [a[i] / b[i] for i in range(len(a))]
mean_t = 0
sigma_t = 1
mean_s = np.mean(r)
sigma_s = np.std(r)
num_bins = 45
n, bins, patches = plt.hist(
    r, num_bins, histtype='step', color='black', density=True)
y_t = ((1 / (np.sqrt(2 * np.pi) * sigma_t)) *
       np.exp(-0.5 * (1 / sigma_t * (bins - mean_t))**2))
y_s = ((1 / (np.sqrt(2 * np.pi) * sigma_s)) *
       np.exp(-0.5 * (1 / sigma_s * (bins - mean_s))**2))
ax.plot(bins, y_t, color='blue')
ax.plot(bins, y_s, color='red')
ax.set_xlabel(
    '(epse_f - epse_s)/|err_pos_epse_f - err_neg_epse_f|')
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
plt.xlim(-0.6, 0.6)
plt.savefig('epse_discrepancy.pdf')
'''


# histogram + scatter plots
def scatter_hist(x, y, ax, ax_histx, ax_histy):

    ax_histx.tick_params(axis="x", labelbottom=True)
    ax_histy.tick_params(axis="y", labelleft=True)

    # the scatter plot:
    ax.scatter(x, y, color='black')
    ax.set_xlabel('n_s')
    ax.set_ylabel('n_f')

    # ax.set_xscale('log')
    ax.set_yscale('log')

    bins_x = 35
    bins_y = 80

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

# scatter_hist(E_s_pow, E_f_pow, ax, ax_histx, ax_histy)

# scatter_hist(Eta0_s, Eta0_f, ax, ax_histx, ax_histy)

# scatter_hist(GammaB_s, GammaB_f, ax, ax_histx, ax_histy)

# scatter_hist(p_s, p_f, ax, ax_histx, ax_histy)

scatter_hist(n_s, n_f_pow, ax, ax_histx, ax_histy)
plt.savefig('n_scatter_dist.png')

# scatter_hist(theta_obs_s, theta_obs_f, ax, ax_histx, ax_histy)

# scatter_hist(epsB_s_pow, epsB_f_pow, ax, ax_histx, ax_histy)

# scatter_hist(epse_s_pow, epse_f_pow, ax, ax_histx, ax_histy)
'''
plt.show()
