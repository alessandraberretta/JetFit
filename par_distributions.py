import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import streamlit as st


E = []
Eta0 = []
GammaB = []
epsB = []
epse = []
n = []
p = []
theta_obs = []
chi2_red = []
Class_T90 = []

Info = {
    # Fitting parameters + chi2 reduced
    'Fit': np.array(['Eta0', 'GammaB', 'theta_obs', 'E', 'epsb', 'epse', 'n', 'p', 'chi2_red', 'Class_T90']),
    # Set parameters in log scale
    'Log': np.array(['E', 'n', 'epse', 'epsb']),
    'LogType': 'Log10'
}


# function that reads csv file
@st.cache
def read_data_file(file, sep='\t'):
    return pd.read_csv(file, sep)


# function for histograms
def hist_pars(par, Info, choosen_pars, s_GRB, l_GRB, T90_classification, num_bins):
    fig, ax = plt.subplots(clear=True)
    if choosen_pars in Info['Log']:
        if T90_classification:
            plt.hist(pow(10, np.array(s_GRB)), bins=num_bins,
                     histtype='step', color='blue')
            plt.hist(pow(10, np.array(l_GRB)), bins=num_bins,
                     histtype='step', color='red')
            mean_s = np.mean(s_GRB)
            sigma_s = np.std(s_GRB)
            mean_l = np.mean(l_GRB)
            sigma_l = np.std(l_GRB)
            textstr_s = '\n'.join((
                r'$\mathrm{mean}=%.2f$' % (mean_s, ),
                r'$\sigma=%.2f$' % (sigma_s, )))
            textstr_l = '\n'.join((
                r'$\mathrm{mean}=%.2f$' % (mean_l, ),
                r'$\sigma=%.2f$' % (sigma_l, )))
            ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='blue')
            ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='red')
        else:
            plt.hist(pow(10, par), bins=num_bins,
                     histtype='step', color='black')
            mean = np.mean(par)
            sigma = np.std(par)
            textstr_s = '\n'.join((
                r'$\mathrm{mean}=%.2f$' % (mean, ),
                r'$\sigma=%.2f$' % (sigma, )))
            ax.text(0.76, 0.99, textstr_s, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='black')
    else:
        if T90_classification:
            plt.hist(s_GRB, bins=num_bins,
                     histtype='step', color='blue')
            plt.hist(l_GRB, bins=num_bins,
                     histtype='step', color='red')
            mean_s = np.mean(s_GRB)
            sigma_s = np.std(s_GRB)
            mean_l = np.mean(l_GRB)
            sigma_l = np.std(l_GRB)
            textstr_s = '\n'.join((
                r'$\mathrm{mean}=%.2f$' % (mean_s, ),
                r'$\sigma=%.2f$' % (sigma_s, )))
            textstr_l = '\n'.join((
                r'$\mathrm{mean}=%.2f$' % (mean_l, ),
                r'$\sigma=%.2f$' % (sigma_l, )))
            ax.text(0.76, 0.99, textstr_s, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='blue')
            ax.text(0.76, 0.89, textstr_l, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='red')
        else:
            plt.hist(par, bins=num_bins, histtype='step', color='black')
            mean = np.mean(par)
            sigma = np.std(par)
            textstr_s = '\n'.join((
                r'$\mathrm{mean}=%.2f$' % (mean, ),
                r'$\sigma=%.2f$' % (sigma, )))
            ax.text(0.76, 0.99, textstr_s, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='black')
    return fig


def main():

    # config streamlit
    st.set_page_config(layout="wide")
    st.write("""
    # Fit parameters distributions
    """)

    choosen_pars = st.sidebar.selectbox(
        'Choose a parameter', ('Eta0', 'GammaB', 'theta_obs', 'E', 'epsb', 'epse', 'n', 'p', 'chi2_red'))
    T90_classification = st.sidebar.checkbox(
        'Long and short classification', True)

    # read data
    path_dirs = '/Users/alessandraberretta/'
    dirs_list = [path_dirs +
                 dir for dir in os.listdir(path_dirs) if dir.endswith('results')]
    summary_list = []
    for elm in dirs_list:
        for file in [elm + '/' + file for file in os.listdir(elm) if file.startswith('new_sum_')]:
            summary_list.append(file)

    for idx, elm in enumerate(summary_list):

        DF = read_data_file(elm)
        E.append(float(DF['Values'][3]))
        Eta0.append(float(DF['Values'][0]))
        GammaB.append(float(DF['Values'][1]))
        epsB.append(float(DF['Values'][4]))
        epse.append(float(DF['Values'][5]))
        n.append(float(DF['Values'][6]))
        p.append(float(DF['Values'][7]))
        theta_obs.append(float(DF['Values'][2]))
        chi2_red.append(float(DF['Values'][12]))
        Class_T90.append(DF['Values'][13])

    Pars = {
        'E': np.array(E),
        'n': np.array(n),
        'Eta0': np.array(Eta0),
        'GammaB': np.array(GammaB),
        'theta_obs': np.array(theta_obs),
        'epse': np.array(epse),
        'epsb': np.array(epsB),
        'p': np.array(p),
        'chi2_red': np.array(chi2_red),
        'Class_T90': np.array(Class_T90)
    }

    inf_in_chi2red = np.isinf(Pars['chi2_red'])

    for idx, elm in enumerate(inf_in_chi2red):
        if elm == True:
            Pars['chi2_red'][idx] = -2

    l_GRB = []
    s_GRB = []

    for idx, elm in enumerate(Pars['Class_T90']):
        if elm == 'S':
            s_GRB.append(Pars[choosen_pars][idx])
        else:
            l_GRB.append(Pars[choosen_pars][idx])

    Latex = {
        'E': r'$log_{10}E_{j,50}$',
        'n': r'$log_{10}n_{0,0}$',
        'Eta0': r'$\eta_0$',
        'GammaB': r'$\gamma_B$',
        'theta_obs': r'$\theta_{obs}$',
        'epse': r'$log_{10}\epsilon_e$',
        'epsb': r'$log_{10}\epsilon_B$',
        'p': r'$p$',
        'chi2_red': r'$\chi^{2}_{red}$',
        'Class_T90': r'$T_{90}$'
    }

    Label = []
    for x in Info['Fit']:
        Label.append(Latex[x])

    hist = hist_pars(Pars[choosen_pars], Info,
                     choosen_pars, s_GRB, l_GRB, T90_classification, num_bins=15)

    st.pyplot(hist)


if __name__ == "__main__":
    main()
