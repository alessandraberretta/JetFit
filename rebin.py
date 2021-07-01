import sys
from tokenize import group
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from decimal import Decimal
import os
import streamlit as st
import statistics
from group_analysis import stdev_groups
from rebin_analysis import rebin_data
from fit_std_analysis import remove_flare
from plot_analysis import lc_plot, scatter_flare


# function that reads csv file
@st.cache
def read_data_file(file):
    return pd.read_csv(file)


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
    group = st.sidebar.checkbox("group the slopes", True)
    y_min = st.sidebar.slider('y_min scatter plot', min_value=float(-500),
                              max_value=float(0), value=float(-100))
    y_max = st.sidebar.slider('y_max scatter plot', min_value=float(0),
                              max_value=float(500), value=float(100))
    sigma_fit_slopes = st.sidebar.number_input('Sigma of slopes distribution', min_value=float(0),
                                               max_value=float(1000), value=float(6.6))
    t_0 = st.sidebar.number_input('Initial cut time', min_value=float(0),
                                  max_value=float(50000), value=float(5000))
    t_rebin = st.sidebar.number_input('Rebin time', min_value=float(
        0), max_value=float(10000), value=float(500))
    group_points = st.sidebar.number_input('Number of points for each group', min_value=int(
        2), max_value=int(1000), value=int(3))

    # read data file
    DF = read_data_file(path_single_GRB)

    lc_original = lc_plot(GRB, DF['Times'].values, DF['Fluxes'].values, DF['FluxErrs'].values, 'black',
                          original=True, rebin=False, flare=False, group=False)
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
                            'green', original=False, rebin=False, flare=True, group=False)
            st.pyplot(scatter)
            col1, col2 = st.beta_columns(2)
            with col1:
                st.pyplot(lc_rebin)
            with col2:
                st.pyplot(lc_nf)
        else:
            if group:
                stdev, mean_stdev_grouped, std_grouped_slopes_removed, lc_group = stdev_groups(
                    GRB, times_mean, fluxes_mean, fluxerrs_mean, group_points)
                Times_red = times_mean
                Fluxes_red = fluxes_mean
                FluxErrs_red = fluxerrs_mean
                col1, col2 = st.beta_columns(2)
                with col1:
                    st.pyplot(lc_original)
                with col2:
                    st.pyplot(lc_rebin)

                st.text("List of std relative to each group")
                st.dataframe(data=stdev, width=None, height=None)

                st.pyplot(lc_group)
    else:
        if flare:
            times_rebin, slopes, times_red, slopes_red_2, Times_red, Fluxes_red, FluxErrs_red = remove_flare(DF['Times'].values, DF['Fluxes'].values,
                                                                                                             DF['FluxErrs'].values, GRB, sigma_fit_slopes)
            scatter = scatter_flare(times_rebin, slopes, times_red,
                                    slopes_red_2, 'black', 'red', y_min, y_max, GRB)
            lc_nf = lc_plot(GRB, Times_red, Fluxes_red, FluxErrs_red,
                            'green', original=False, rebin=False, flare=True, group=False)
            st.pyplot(scatter)
            col1, col2 = st.beta_columns(2)
            with col1:
                st.pyplot(lc_original)
            with col2:
                st.pyplot(lc_nf)
        else:
            Times_red = DF['Times'].values
            Fluxes_red = DF['Fluxes'].values
            FluxErrs_red = DF['FluxErrs'].values

    d2 = {'Times': Times_red, 'Fluxes': Fluxes_red, 'FluxErrs': FluxErrs_red}
    df2 = pd.DataFrame(data=d2, columns=['Times', 'Fluxes', 'FluxErrs'])
    df2["Times"] = ['%.6E' % Decimal(x) for x in df2['Times']]
    df2["Fluxes"] = ['%.6E' % Decimal(y) for y in df2['Fluxes']]
    df2["FluxErrs"] = ['%.6E' % Decimal(z) for z in df2['FluxErrs']]
    df2.to_csv(GRB + "_nf" + ".csv",  sep='\t')

    path_lc = '/Users/alessandraberretta/JetFit/2013_rebin_removeflare_results/'
    list_fitted_GRB = [path_lc +
                       file for file in os.listdir(path_lc) if file.startswith('lc_')]
    list_summary_fitted_GRB = [path_lc +
                               file for file in os.listdir(path_lc) if file.startswith('summary_')]
    chi2red_list = 0

    for summary in list_summary_fitted_GRB:
        if GRB in summary.split('/')[5]:
            df = pd.read_csv(summary, sep='\t')
            chi2red_list = df['Values'][12]
    for fittedGRB in list_fitted_GRB:
        if GRB in fittedGRB:
            st.image(
                fittedGRB, caption=GRB + " " + " " + " " + "chi2_red = " +
                str(chi2red_list))


if __name__ == "__main__":
    main()
