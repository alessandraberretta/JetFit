import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from decimal import Decimal
import os


# read data
path_GRB = '/Users/alessandraberretta/JetFit/2016/2016def/'
list_file = [path_GRB +
             file for file in os.listdir(path_GRB) if file.startswith('GRB_')]

single = True
path_single_GRB = '/Users/alessandraberretta/JetFit/2013/2013def/GRB_130907A_1.238_def.csv'
GRB = path_single_GRB[path_single_GRB.rfind(
    '/')+1:path_single_GRB.rfind('_def')]

path_out_nolog10 = '/Users/alessandraberretta/JetFit/flare_analysis/2016/no_log10/'
path_out_log10 = '/Users/alessandraberretta/JetFit/flare_analysis/2016/log10/'

path_out_nolog10_half = '/Users/alessandraberretta/JetFit/flare_analysis/2016/no_log10/half/'
path_out_log10_half = '/Users/alessandraberretta/JetFit/flare_analysis/2016/log10/half/'


def scatter_flare(x, y, color, ymin, ymax, flare=False):
    plt.clf()
    fig, ax = plt.subplots()
    plt.scatter(x, y, c=color)
    ax.set_xscale('log')
    ax.set_ylim(ymin, ymax)
    if single:
        if flare:
            plt.savefig("slopes_times_" + GRB + ".pdf")
        else:
            plt.savefig("slopes_times_" + GRB + "_nf.pdf")
    else:
        if flare:
            plt.savefig(path_out_nolog10_half + "slopes_times_" +
                        elm[elm.rfind('/') + 1:elm.rfind('_def')] + ".pdf")
        else:
            plt.savefig(path_out_nolog10_half + "slopes_times_" +
                        elm[elm.rfind('/') + 1:elm.rfind('_def')] + "_nf.pdf")


def lc_no_flare(times, fluxes, fluxes_err, color):
    plt.clf()
    fig, ax = plt.subplots()
    ax.errorbar(times, fluxes, yerr=fluxes_err, fmt='.', color=color)
    ax.set_yscale('log')
    ax.set_xscale('log')
    if single:
        plt.savefig("lc_" + GRB + "_nf.pdf")
    else:
        plt.savefig(path_out_nolog10_half + "lc_" +
                    elm[elm.rfind('/') + 1:elm.rfind('_def')] + "_nf.pdf")


# sigma_dist_slopes = 29/2
# sigma_dist_slopes = (2.312e-13)/2
# sigma_singleGRB_slopes = 9.394e-13/2
sigma_singleGRB_slopes = 100

slopes = []
intercepts = []

if single:

    DF = pd.read_csv(path_single_GRB)
    Times, Fluxes, FluxErrs = DF['Times'].values, DF['Fluxes'].values, DF['FluxErrs'].values

    for idx in range(len(Times)-1):
        slope, intercept, _, _, _ = linregress([np.log10(Times[idx]), np.log10(
            Times[idx+1])], [np.log10(Fluxes[idx]), np.log10(Fluxes[idx+1])])
        # slope, intercept, _, _, _ = linregress(
        # [Times[idx], Times[idx+1]], [Fluxes[idx], Fluxes[idx+1]])
        slopes.append(slope)
        intercepts.append(intercept)

    # compute the times_mean and the first/last point of each bin to plot times vs slopes
    Times_mean = []
    for idx in range(len(Times)-1):
        Times_mean.append(Times[idx])
    first_point_t = []
    last_point_t = []
    first_point_f = []
    last_point_f = []
    for idx in range(len(Times)-1):
        first_point_t.append(Times[idx])
        last_point_t.append(Times[idx+1])
        first_point_f.append(Fluxes[idx])
        last_point_f.append(Fluxes[idx+1])

    d = {'Slopes': slopes, 'Intercepts': intercepts, 'Times_mean': Times_mean, 'First_time': first_point_t,
         'Last_time': last_point_t, 'First_flux': first_point_f, 'Last_flux': last_point_f}
    df = pd.DataFrame(data=d, columns=[
        'Slopes', 'Intercepts', 'Times_mean', 'First_time', 'Last_time', 'First_flux', 'Last_flux'])
    df["Slopes"] = ['%.6E' % Decimal(x) for x in df['Slopes']]
    df["Intercepts"] = ['%.6E' % Decimal(y) for y in df['Intercepts']]
    df["Times_mean"] = ['%.6E' % Decimal(z) for z in df['Times_mean']]
    df["First_time"] = ['%.6E' % Decimal(z) for z in df['First_time']]
    df["Last_time"] = ['%.6E' % Decimal(z) for z in df['Last_time']]
    df["First_flux"] = ['%.6E' % Decimal(z) for z in df['First_flux']]
    df["Last_flux"] = ['%.6E' % Decimal(z) for z in df['Last_flux']]
    df.to_csv("info_" + GRB + ".csv",  sep='\t')

    scatter_flare(Times_mean, slopes, 'black', -100, 100)
    # scatter_flare(Times_mean, slopes, 'black', -20e-12, 20e-12)

    slopes_red = []
    dropped_idx = []
    for idx, val in enumerate(slopes):
        print(val, sigma_singleGRB_slopes)
        # if val > sigma_singleGRB_slopes or val < -sigma_singleGRB_slopes:
        if abs(val) > sigma_singleGRB_slopes:
            print('no flare')
            slopes_red.append(val)
            i = df.loc[df['Slopes'] == '%.6E' % Decimal(val)].index
            '''
            if val > 0:
                i = df.loc[i+1].index
            '''
            dropped_idx.append(i[0])
    dropped_idx_arr = np.unique(np.asarray(dropped_idx))
    for k in dropped_idx_arr:
        df = df.drop(k)
    times_mean_red = [float(i) for i in df['Times_mean'].values]
    slopes_red_2 = [float(i) for i in df['Slopes'].values]

    scatter_flare(times_mean_red, slopes_red_2, 'red', -100, 100, flare=True)
    # scatter_flare(times_mean_red, slopes_red_2,
    # 'red', -20e-12, 20e-12, flare=True)

    Times_red = np.delete(Times, dropped_idx)
    Fluxes_red = np.delete(Fluxes, dropped_idx)
    FluxErrs_red = np.delete(FluxErrs, dropped_idx)

    lc_no_flare(Times_red, Fluxes_red, FluxErrs_red, 'red')

else:

    for elm in list_file:

        DF = pd.read_csv(elm)
        Times, Fluxes, FluxErrs = DF['Times'].values, DF['Fluxes'].values, DF['FluxErrs'].values

        # get slopes and intercepts of binned data
        slopes = []
        intercepts = []

        for idx in range(len(Times)-1):
            # slope, intercept, _, _, _ = linregress([np.log10(Times[idx]), np.log10(
            # Times[idx+1])], [np.log10(Fluxes[idx]), np.log10(Fluxes[idx+1])])
            slope, intercept, _, _, _ = linregress(
                [Times[idx], Times[idx+1]], [Fluxes[idx], Fluxes[idx+1]])
            slopes.append(slope)
            intercepts.append(intercept)

        # compute the times_mean and the first/last point of each bin to plot times vs slopes
        Times_mean = []
        for idx in range(len(Times)-1):
            Times_mean.append(Times[idx])
        first_point_t = []
        last_point_t = []
        first_point_f = []
        last_point_f = []
        for idx in range(len(Times)-1):
            first_point_t.append(Times[idx])
            last_point_t.append(Times[idx+1])
            first_point_f.append(Fluxes[idx])
            last_point_f.append(Fluxes[idx+1])

        d = {'Slopes': slopes, 'Intercepts': intercepts, 'Times_mean': Times_mean, 'First_time': first_point_t,
             'Last_time': last_point_t, 'First_flux': first_point_f, 'Last_flux': last_point_f}
        df = pd.DataFrame(data=d, columns=[
            'Slopes', 'Intercepts', 'Times_mean', 'First_time', 'Last_time', 'First_flux', 'Last_flux'])
        df["Slopes"] = ['%.6E' % Decimal(x) for x in df['Slopes']]
        df["Intercepts"] = ['%.6E' % Decimal(y) for y in df['Intercepts']]
        df["Times_mean"] = ['%.6E' % Decimal(z) for z in df['Times_mean']]
        df["First_time"] = ['%.6E' % Decimal(z) for z in df['First_time']]
        df["Last_time"] = ['%.6E' % Decimal(z) for z in df['Last_time']]
        df["First_flux"] = ['%.6E' % Decimal(z) for z in df['First_flux']]
        df["Last_flux"] = ['%.6E' % Decimal(z) for z in df['Last_flux']]
        df.to_csv("info_" + elm[elm.rfind('/') +
                                1:elm.rfind('_def')] + ".csv",  sep='\t')

        # scatter_flare(Times_mean, slopes, 'black', -100, 100)
        scatter_flare(Times_mean, slopes, 'black', -20e-12, 20e-12)

        slopes_red = []
        dropped_idx = []
        for idx, val in enumerate(slopes):
            if val > sigma_dist_slopes or val < -sigma_dist_slopes:
                slopes_red.append(val)
                i = df.loc[df['Slopes'] == '%.6E' % Decimal(val)].index
                if val > 0:
                    i = df.loc[i+1].index
                dropped_idx.append(i[0])
        dropped_idx_arr = np.unique(np.asarray(dropped_idx))
        for k in dropped_idx_arr:
            df = df.drop(k)
        times_mean_red = [float(i) for i in df['Times_mean'].values]
        slopes_red_2 = [float(i) for i in df['Slopes'].values]

        # scatter_flare(times_mean_red, slopes_red_2, 'red', -100, 100, flare=True)
        scatter_flare(times_mean_red, slopes_red_2,
                      'red', -20e-12, 20e-12, flare=True)

        Times_red = np.delete(Times, dropped_idx)
        Fluxes_red = np.delete(Fluxes, dropped_idx)
        FluxErrs_red = np.delete(FluxErrs, dropped_idx)

        lc_no_flare(Times_red, Fluxes_red, FluxErrs_red, 'red')


'''
GRB_slopes = 'slopes_info_131103A_0.599.csv'
DF = pd.read_csv(GRB_slopes, sep='\t')
slopes_GRB = DF['Slopes'].values


# select and drop the slopes that are outside 1-sigma distribution of year slopes
slopes_red = []
dropped_idx = []
print(np.std(slopes_GRB))
for idx, val in enumerate(slopes):
    # if val > (1.40e-12)/2 or val < (-1.40e-12)/2:
    if val > 29 or val < -29:
        slopes_red.append(val)
        i = df.loc[df['Slopes'] == '%.6E' % Decimal(val)].index
        #print(f'Actual index: {i[0]} - {val}')
        if val > 0:
            i = df.loc[i+1].index
        dropped_idx.append(i[0])
dropped_idx_arr = np.unique(np.asarray(dropped_idx))
print(dropped_idx_arr)
for k in dropped_idx_arr:
    df = df.drop(k)


# plot times vs slopes of reduced data and write a csv file with reduced info
times_mean_red = [float(i) for i in df['Times_mean'].values]
slopes_red_2 = [float(i) for i in df['Slopes'].values]
plt.scatter(times_mean_red, slopes_red_2, c='red')
# ax.set_ylim(-25, 25)
ax.set_xscale('log')
plt.savefig("slopes_times_mean_131103A_0.599_nf.pdf")
# plt.show()


# drop the corresponding times and fluxes to get the reduced lc with "no flares"
Times_red = np.delete(Times, dropped_idx)
Fluxes_red = np.delete(Fluxes, dropped_idx)
FluxErrs_red = np.delete(FluxErrs, dropped_idx)
ax.errorbar(Times_red, Fluxes_red, yerr=FluxErrs_red, fmt='.', color='red')
ax.set_yscale('log')
ax.set_xscale('log')
plt.savefig('lc_130420A_131103A_0.599.pdf')
'''
