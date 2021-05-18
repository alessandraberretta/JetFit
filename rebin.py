import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from decimal import Decimal
import os
import streamlit as st
import base64


st.set_page_config(layout="wide")
st.write("""
# Rebin analysis on a sample of *Swift*-XRT GRBs
""")

st.sidebar.write("Analysis flags")
single = st.sidebar.checkbox("single", True)
# flare = st.sidebar.checkbox("flare", True)
partial_rebin = st.sidebar.checkbox("partial_rebin", True)
original_lc = st.sidebar.checkbox("original_lc", True)

# flag for removed flare
flare = True

# flag for single GRB
# single = True
# flag for removed flare
# flare = False
# rebin only in the final points of the lc
# partial_rebin = True


# read data
GRB_list = ['/Users/alessandraberretta/JetFit/2013/2013def/GRB_130420A_1.297_def.csv',
            '/Users/alessandraberretta/JetFit/2013/2013def/GRB_130505A_2.27_def.csv',
            '/Users/alessandraberretta/JetFit/2013/2013def/GRB_130606A_5.91_def.csv',
            '/Users/alessandraberretta/JetFit/2013/2013def/GRB_130907A_1.238_def.csv',
            '/Users/alessandraberretta/JetFit/2014/2014def/GRB_140304A_5.28_def.csv',
            '/Users/alessandraberretta/JetFit/2015/2015def/GRB_150206A_2.087_def.csv',
            '/Users/alessandraberretta/JetFit/2015/2015def/GRB_150821A_0.755_def.csv',
            '/Users/alessandraberretta/JetFit/2016/2016def/GRB_160117B_0.87_def.csv']
choosen_GRB = st.sidebar.selectbox(
    'GRB file', ('GRB 130420A', 'GRB 130505A', 'GRB 130606A', 'GRB 130907A', 'GRB 140304A', 'GRB 150206A',
                 'GRB 150821A', 'GRB 160117B'))
for file in GRB_list:
    if choosen_GRB[choosen_GRB.rfind(' ')+1:] in file:
        path_single_GRB = file
# path_single_GRB = st.sidebar.text_input('GRB file', '/Users/alessandraberretta/JetFit/2013/2013def/GRB_130420A_1.297_def.csv')
# path_single_GRB = '/Users/alessandraberretta/JetFit/2015/2015def/GRB_150821A_0.755_def.csv'
GRB = path_single_GRB[path_single_GRB.rfind(
    '/')+1:path_single_GRB.rfind('_def')]


@st.cache
def read_data_file(file):
    return pd.read_csv(file)


def scatter_flare(x_1, y_1, x_2, y_2, color_1, color_2, ymin, ymax):
    fig, ax = plt.subplots(clear=True)
    # plt.clf()
    plt.scatter(x_1, y_1, c=color_1)
    plt.scatter(x_2, y_2, c=color_2)
    ax.set_xscale('log')
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('times_mean (s)')
    ax.set_ylabel('slopes')
    if single:
        plt.title('Slopes vs Times_mean scatter plot')
        # plt.savefig("slopes_times_" + GRB + "_rebin.pdf")
    else:
        plt.savefig("slopes_times_" + GRB + "_nf.pdf")
    return fig


def lc_plot(times, fluxes, fluxes_err, color, original_lc=True, flare=True):
    fig, ax = plt.subplots(clear=True)
    # plt.clf()
    ax.errorbar(times, fluxes, yerr=fluxes_err, fmt='.', color=color)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('times (s)')
    ax.set_ylabel('flux (erg/cm^2/s)')
    if single:
        if original_lc:
            plt.title('Original lightcurve')
            # plt.savefig("lc_" + GRB + ".pdf")
        if flare:
            plt.title('Re-binned lightcurve')
            # plt.savefig("lc_" + GRB + "_rebin.pdf")
        if original_lc == False and flare == False:
            plt.title('Removed flares lightcurve')
            # plt.savefig("lc_" + GRB + "_nf.pdf")
    else:
        plt.savefig("lc_" + GRB + "_nf.pdf")
    return fig


def weighted_error_mean(var, err_var):
    return np.average(var, axis=0, weights=[1/i for i in err_var])


temp_times = []
temp_fluxes = []
temp_fluxerrs = []

times_mean = []
fluxes_mean = []
fluxerrs_mean = []

t_0 = st.sidebar.slider('Initial cut time', min_value=float("{:e}".format(0)),
                        max_value=float("{:e}".format(50000)), value=float("{:e}".format(5000)))
# t_0 = 5e3

sigma_fit_slopes = st.sidebar.slider('Sigma of slopes distribution', min_value=float("{:e}".format(1)),
                                     max_value=float("{:e}".format(50)), value=float("{:e}".format(6.6)))
# sigma_fit_slopes = 10

t_rebin = st.sidebar.slider('Rebin time', min_value=float(0),
                            max_value=float(2000), value=float(500))
# t_rebin = 300

if single:

    DF = read_data_file(path_single_GRB)
    Times, Fluxes, FluxErrs = DF['Times'].values, DF['Fluxes'].values, DF['FluxErrs'].values
    Times_0 = DF['Times'][0]

    for idx, val in enumerate(Times):
        if val > t_0:
            if len(temp_times) == 0:
                temp_times.append(val)
                temp_fluxes.append(Fluxes[idx])
                temp_fluxerrs.append(FluxErrs[idx])
            else:
                if (val - temp_times[0]) < t_rebin:
                    temp_times.append(val)
                    temp_fluxes.append(Fluxes[idx])
                    temp_fluxerrs.append(FluxErrs[idx])
                else:
                    # print(val - temp_times[0], val, Fluxes[idx])
                    times_mean.append(np.mean(temp_times))
                    fluxes_mean.append(weighted_error_mean(
                        temp_fluxes, temp_fluxerrs))
                    fluxerrs_mean.append(weighted_error_mean(
                        temp_fluxerrs, temp_fluxerrs))
                    temp_times.clear()
                    temp_fluxes.clear()
                    temp_fluxerrs.clear()
                    temp_times.append(val)
                    temp_fluxes.append(Fluxes[idx])
                    temp_fluxerrs.append(FluxErrs[idx])
        elif val >= Times_0 and val < t_0:
            if partial_rebin:
                times_mean.append(val)
                fluxes_mean.append(Fluxes[idx])
                fluxerrs_mean.append(FluxErrs[idx])

    lc_original = lc_plot(Times, Fluxes, FluxErrs, 'black',
                          original_lc=True, flare=False)

    lc_rebin = lc_plot(times_mean, fluxes_mean,
                       fluxerrs_mean, 'red', original_lc=False)

    slopes = []
    intercepts = []
    Times_mean_rebin = []
    first_point_t = []
    last_point_t = []
    first_point_f = []
    last_point_f = []

    for idx in range(len(times_mean)-1):
        # slope, intercept, _, _, _ = linregress([np.log10(Times[idx]), np.log10(
        # Times[idx+1])], [np.log10(Fluxes[idx]), np.log10(Fluxes[idx+1])])
        slope, intercept, _, _, _ = linregress([np.log10(times_mean[idx]), np.log10(
            times_mean[idx+1])], [np.log10(fluxes_mean[idx]), np.log10(fluxes_mean[idx+1])])
        slopes.append(slope)
        intercepts.append(intercept)
        Times_mean_rebin.append(times_mean[idx])
        first_point_t.append(times_mean[idx])
        last_point_t.append(times_mean[idx+1])
        first_point_f.append(fluxes_mean[idx])
        last_point_f.append(fluxes_mean[idx+1])

    d = {'Slopes': slopes, 'Intercepts': intercepts, 'Times_mean': Times_mean_rebin, 'First_time': first_point_t,
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

    slopes_red = []
    dropped_idx = []
    for idx, val in enumerate(slopes):
        # if val > sigma_singleGRB_slopes or val < -sigma_singleGRB_slopes:
        if abs(val) > sigma_fit_slopes:
            slopes_red.append(val)
            i = df.loc[df['Slopes'] == '%.6E' % Decimal(val)].index
            dropped_idx.append(i[0])
    dropped_idx_arr = np.unique(np.asarray(dropped_idx))
    for k in dropped_idx_arr:
        df = df.drop(k)
    times_mean_red = [float(i) for i in df['Times_mean'].values]
    slopes_red_2 = [float(i) for i in df['Slopes'].values]

    scatter = scatter_flare(Times_mean_rebin, slopes, times_mean_red,
                            slopes_red_2, 'black', 'red', -200, 200)

    st.pyplot(scatter)
    st.pyplot(lc_original)

    Times_red = np.delete(times_mean, dropped_idx)
    Fluxes_red = np.delete(fluxes_mean, dropped_idx)
    FluxErrs_red = np.delete(fluxerrs_mean, dropped_idx)

    lc_nf = lc_plot(Times_red, Fluxes_red, FluxErrs_red,
                    'green', original_lc=False, flare=False)

    col1, col2 = st.beta_columns(2)
    if original_lc:
        with col1:
            st.pyplot(lc_rebin)
        with col2:
            st.pyplot(lc_nf)
    else:
        st.pyplot(lc_nf)

    # st.balloons()