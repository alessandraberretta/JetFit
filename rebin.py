import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from decimal import Decimal
import os
import streamlit as st


# function that reads csv file
@st.cache
def read_data_file(file):
    return pd.read_csv(file)


# function for scatter plot
def scatter_flare(x_1, y_1, x_2, y_2, color_1, color_2, ymin, ymax, GRB):
    fig, ax = plt.subplots(clear=True)
    # plt.clf()
    plt.scatter(x_1, y_1, c=color_1)
    plt.scatter(x_2, y_2, c=color_2)
    ax.set_xscale('log')
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel('times (s)')
    ax.set_ylabel('slopes')
    plt.title('Slopes vs times scatter plot')
    plt.savefig("slopes_times_" + GRB + "_rebin.pdf")
    return fig


# function for lighcurve plot
def lc_plot(GRB, times, fluxes, fluxes_err, color, original=True, rebin=True, flare=True):
    fig, ax = plt.subplots(clear=True)
    # plt.clf()
    ax.errorbar(times, fluxes, yerr=fluxes_err, fmt='.', color=color)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('times (s)')
    ax.set_ylabel('flux (erg/cm^2/s)')
    if original:
        plt.title('Original lightcurve')
        plt.savefig("lc_" + GRB + ".pdf")
    if rebin:
        plt.title('Re-binned lightcurve')
        plt.savefig("lc_" + GRB + "_rebin.pdf")
    if flare:
        plt.title('Removed flares lightcurve')
        plt.savefig("lc_" + GRB + "_nf.pdf")
    return fig


# function for weighted mean computation
def weighted_error_mean(var, err_var):
    return np.average(var, axis=0, weights=[1/i for i in err_var])


# function for data rebinning
def rebin_data(GRB, DF, t_0, t_rebin, partial_rebin):
    # define temporary/empty arrays for the following analysis
    temp_times = []
    temp_fluxes = []
    temp_fluxerrs = []
    times_mean = []
    fluxes_mean = []
    fluxerrs_mean = []

    Times, Fluxes, FluxErrs = DF['Times'].values, DF['Fluxes'].values, DF['FluxErrs'].values
    Times_0 = DF['Times'][0]

    for idx, val in enumerate(Times):
        if val > t_0:
            if len(temp_times) == 0:
                temp_times.append(val)
                temp_fluxes.append(Fluxes[idx])
                temp_fluxerrs.append(FluxErrs[idx])
            else:
                if abs(val - temp_times[0]) < t_rebin:
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

    # get the rebinned lightcurve
    lc_rebin = lc_plot(GRB, times_mean, fluxes_mean,
                       fluxerrs_mean, 'red', original=False, rebin=True, flare=False)

    return lc_rebin, times_mean, fluxes_mean, fluxerrs_mean


# function for remove flares from data
def remove_flare(times, fluxes, fluxerrs, GRB, sigma_fit_slopes):

    slopes = []
    intercepts = []
    times_rebin = []
    first_point_t = []
    last_point_t = []
    first_point_f = []
    last_point_f = []

    # compute the slopes between couples of points
    for idx in range(len(times)-1):
        slope, intercept, _, _, _ = linregress([np.log10(times[idx]), np.log10(
            times[idx+1])], [np.log10(fluxes[idx]), np.log10(fluxes[idx+1])])
        slopes.append(slope)
        intercepts.append(intercept)
        times_rebin.append(times[idx])
        first_point_t.append(times[idx])
        last_point_t.append(times[idx+1])
        first_point_f.append(fluxes[idx])
        last_point_f.append(fluxes[idx+1])

    # create a dataframe with all the information
    d = {'Slopes': slopes, 'Intercepts': intercepts, 'times': times_rebin, 'First_time': first_point_t,
         'Last_time': last_point_t, 'First_flux': first_point_f, 'Last_flux': last_point_f}
    df = pd.DataFrame(data=d, columns=[
        'Slopes', 'Intercepts', 'times', 'First_time', 'Last_time', 'First_flux', 'Last_flux'])
    df["Slopes"] = ['%.6E' % Decimal(x) for x in df['Slopes']]
    df["Intercepts"] = ['%.6E' % Decimal(y) for y in df['Intercepts']]
    df["times"] = ['%.6E' % Decimal(z) for z in df['times']]
    df["First_time"] = ['%.6E' % Decimal(z) for z in df['First_time']]
    df["Last_time"] = ['%.6E' % Decimal(z) for z in df['Last_time']]
    df["First_flux"] = ['%.6E' % Decimal(z) for z in df['First_flux']]
    df["Last_flux"] = ['%.6E' % Decimal(z) for z in df['Last_flux']]
    df.to_csv("info_" + GRB + ".csv",  sep='\t')

    # using the sigma of the slopes distribution to filter the data
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
    times_red = [float(i) for i in df['times'].values]
    slopes_red_2 = [float(i) for i in df['Slopes'].values]

    # delete from the original arrays the data points that not survived the sigma cut
    Times_red = np.delete(times, dropped_idx)
    Fluxes_red = np.delete(fluxes, dropped_idx)
    FluxErrs_red = np.delete(fluxerrs, dropped_idx)

    return times_rebin, slopes, times_red, slopes_red_2, Times_red, Fluxes_red, FluxErrs_red


def main():

    # config streamlit
    st.set_page_config(layout="wide")
    st.write("""
    # Rebin analysis on a sample of *Swift*-XRT GRBs
    """)

    # read data
    path_GRB = '/Users/alessandraberretta/JetFit/2013/2013def/'
    GRB_list = [path_GRB +
                file for file in os.listdir(path_GRB) if file.startswith('GRB_')]

    choosen_GRB = st.sidebar.selectbox(
        'GRB file', ('GRB 130408A', 'GRB 130418A', 'GRB 130420A', 'GRB 130427B', 'GRB 130505A', 'GRB 130511A',
                     'GRB 130603B', 'GRB 130606A', 'GRB 130610A', 'GRB 130612A', 'GRB 130701A', 'GRB 130702A',
                     'GRB 130831A', 'GRB 130907A', 'GRB 130925A', 'GRB 131004A', 'GRB 131030A', 'GRB 131103A',
                     'GRB 131105A', 'GRB 131108A', 'GRB 131117A', 'GRB 131231A'))
    for file in GRB_list:
        if choosen_GRB[choosen_GRB.rfind(' ')+1:] in file:
            path_single_GRB = file
    GRB = path_single_GRB[path_single_GRB.rfind(
        '/')+1:path_single_GRB.rfind('_def')]

    # some useful sliders and streamlit buttons
    st.sidebar.write("Analysis buttons")
    flare = st.sidebar.checkbox("a flare is present", True)
    partial_rebin = st.sidebar.checkbox("partial_rebin", True)
    # original_lc = st.sidebar.checkbox("original_lc", True)
    y_min = st.sidebar.slider('y_min scatter plot', min_value=float(-250),
                              max_value=float(0), value=float(-100))
    y_max = st.sidebar.slider('y_max scatter plot', min_value=float(0),
                              max_value=float(250), value=float(100))
    sigma_fit_slopes = st.sidebar.number_input('Sigma of slopes distribution', min_value=float(0),
                                               max_value=float(1000), value=float(6.6))
    t_0 = st.sidebar.number_input('Initial cut time', min_value=float(0),
                                  max_value=float(50000), value=float(5000))
    t_rebin = st.sidebar.number_input('Rebin time', min_value=float(
        0), max_value=float(10000), value=float(500))

    # read data file
    DF = read_data_file(path_single_GRB)

    lc_original = lc_plot(GRB, DF['Times'], DF['Fluxes'], DF['FluxErrs'], 'black',
                          original=True, rebin=False, flare=False)
    st.pyplot(lc_original)

    if t_rebin:
        lc_rebin, times_mean, fluxes_mean, fluxerrs_mean = rebin_data(
            GRB, DF, t_0, t_rebin, partial_rebin)
        if flare:
            times_rebin, slopes, times_red, slopes_red_2, Times_red, Fluxes_red, FluxErrs_red = remove_flare(times_mean, fluxes_mean,
                                                                                                             fluxerrs_mean, GRB, sigma_fit_slopes)
            scatter = scatter_flare(times_rebin, slopes, times_red,
                                    slopes_red_2, 'black', 'red', y_min, y_max, GRB)
            lc_nf = lc_plot(GRB, Times_red, Fluxes_red, FluxErrs_red,
                            'green', original=False, rebin=False, flare=True)
            st.pyplot(scatter)
            # plot the two lightcurves with streamlit
            col1, col2 = st.beta_columns(2)
            with col1:
                st.pyplot(lc_rebin)
            with col2:
                st.pyplot(lc_nf)
        else:
            Times_red = times_mean
            Fluxes_red = fluxes_mean
            FluxErrs_red = fluxerrs_mean
            col1, col2 = st.beta_columns(2)
            with col1:
                st.pyplot(lc_original)
            with col2:
                st.pyplot(lc_rebin)
    else:
        if flare:
            times_rebin, slopes, times_red, slopes_red_2, Times_red, Fluxes_red, FluxErrs_red = remove_flare(DF['Times'], DF['Fluxes'],
                                                                                                             DF['FluxErrs'], GRB, sigma_fit_slopes)
            scatter = scatter_flare(times_rebin, slopes, times_red,
                                    slopes_red_2, 'black', 'red', y_min, y_max, GRB)
            lc_nf = lc_plot(GRB, Times_red, Fluxes_red, FluxErrs_red,
                            'green', original=False, rebin=False, flare=True)
            st.pyplot(scatter)
            # plot the two lightcurves with streamlit
            col1, col2 = st.beta_columns(2)
            with col1:
                st.pyplot(lc_original)
            with col2:
                st.pyplot(lc_nf)
        else:
            Times_red = DF['Times']
            Fluxes_red = DF['Fluxes']
            FluxErrs_red = DF['FluxErrs']

    d2 = {'Times': Times_red, 'Fluxes': Fluxes_red, 'FluxErrs': FluxErrs_red}
    df2 = pd.DataFrame(data=d2, columns=['Times', 'Fluxes', 'FluxErrs'])
    df2["Times"] = ['%.6E' % Decimal(x) for x in df2['Times']]
    df2["Fluxes"] = ['%.6E' % Decimal(y) for y in df2['Fluxes']]
    df2["FluxErrs"] = ['%.6E' % Decimal(z) for z in df2['FluxErrs']]
    df2.to_csv(GRB + "_nf" + ".csv",  sep='\t')


if __name__ == "__main__":
    main()
