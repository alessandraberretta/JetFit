import sys
from matplotlib import patches
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib import histograms
import pandas as pd
import os
import streamlit as st
import math


E = []
Eta0 = []
GammaB = []
epsB = []
epse = []
n = []
p = []
theta_obs = []
redshift = []
chi2_red = []
Class_T90 = []
names_dist2 = []


Info = {
    # Fitting parameters + chi2 reduced
    'Fit': np.array(['redshift', 'Eta0', 'GammaB', 'theta_0', 'theta_obs', 'E', 'epsb', 'epse', 'n', 'p', 'thetaobs_theta0', 'GAMMA', 'chi2_red', 'Class_T90']),
    # Set parameters in log scale
    'Log': np.array(['E', 'n', 'epse', 'epsb']),
    'LogType': 'Log10'
}


# function that reads csv file
@st.cache
def read_data_file(file, sep='\t'):
    return pd.read_csv(file, sep)


# function for histograms
def hist_pars(class_by_tSNE, class_by_gammaB, no_classification, class_1_temp, class_2_temp, prob_dist, histo_classification,
              class_small, class_big, par, Info, choosen_pars, theta_obs, theta_obs_s, theta_obs_l, s_GRB, l_GRB, reasonable_long, reasonable_short,
              T90_classification, num_bins, num_bins_class_small, num_bins_class_big):
    # num_bins_s, num_bins_l,
    fig, ax = plt.subplots(clear=True)
    if histo_classification:
        class_1_temp = np.asarray(class_1_temp)
        class_2_temp = np.asarray(class_2_temp)
        if 'theta_obs' or 'thetaobs_theta0' in choosen_pars:
            w = 0.02
            n_1 = math.ceil((class_1_temp.max() - class_1_temp.min())/w)
            n_2 = math.ceil((class_2_temp.max() - class_2_temp.min())/w)
            plt.hist(class_1_temp, bins=n_1,
                     color='red', histtype='step', linewidth=1.5, alpha=0.7, edgecolor='red')
            plt.hist(class_2_temp, bins=n_2,
                     color='blue', histtype='step', linewidth=1.5, alpha=0.7, edgecolor='blue')
            median_1 = np.median(class_1_temp)
            sigma_1 = np.std(class_1_temp)
            median_2 = np.median(class_2_temp)
            sigma_2 = np.std(class_2_temp)
            textstr_1 = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_1, ),
                r'$\sigma=%.2f$' % (sigma_1, )))
            textstr_2 = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_2, ),
                r'$\sigma=%.2f$' % (sigma_2, )))
            ax.text(0.70, 0.99, textstr_1, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            ax.text(0.70, 0.89, textstr_2, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
        else:
            w = 0.3
            n_1 = math.ceil((class_1_temp.max() - class_1_temp.min())/w)
            n_2 = math.ceil((class_2_temp.max() - class_2_temp.min())/w)
            plt.hist(class_1_temp, bins=n_1,
                     color='red', histtype='step', linewidth=1.5, alpha=0.7, edgecolor='red')
            plt.hist(class_2_temp, bins=n_2,
                     color='blue', histtype='step', linewidth=1.5, alpha=0.7, edgecolor='blue')
            median_1 = np.median(class_1_temp)
            sigma_1 = np.std(class_1_temp)
            median_2 = np.median(class_2_temp)
            sigma_2 = np.std(class_2_temp)
            textstr_1 = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_1, ),
                r'$\sigma=%.2f$' % (sigma_1, )))
            textstr_2 = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_2, ),
                r'$\sigma=%.2f$' % (sigma_2, )))
            ax.text(0.70, 0.99, textstr_1, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            ax.text(0.70, 0.89, textstr_2, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
    if prob_dist:
        if 'thetaobs_theta0' in choosen_pars:
            if T90_classification:
                s_GRB = np.asarray(s_GRB)
                l_GRB = np.asarray(l_GRB)
                w = 0.04
                n_s = math.ceil((s_GRB.max() - s_GRB.min())/w)
                n_l = math.ceil((l_GRB.max() - l_GRB.min())/w)
                weight_s = 1/(2*np.pi*np.sin(theta_obs_s))
                weight_l = 1/(2*np.pi*np.sin(theta_obs_l))
                plt.hist(s_GRB, bins=n_s, weights=weight_s,
                         color='blue', histtype='step', linewidth=1.5, alpha=0.7, edgecolor='blue')
                plt.hist(l_GRB, bins=n_l,
                         weights=weight_l, histtype='step', linewidth=1.5, color='red', alpha=0.7, edgecolor='red')
                median_s = np.median(s_GRB)
                sigma_s = np.std(s_GRB)
                median_l = np.median(l_GRB)
                sigma_l = np.std(l_GRB)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_s, ),
                    r'$\sigma=%.2f$' % (sigma_s, )))
                textstr_l = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_l, ),
                    r'$\sigma=%.2f$' % (sigma_l, )))
                ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.6)
                ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes, fontsize=13,
                        verticalalignment='top', color='red', alpha=0.6)
            else:
                weight = 1/(2*np.pi*np.sin(theta_obs))
                n, bins, patches = plt.hist(par, bins=num_bins, color='blue',
                                            weights=weight, alpha=0.3)
                print(n)
                print(bins)
                plt.title(
                    'Probability of observing a GRB with a definite opening angle')
                plt.show()
                median = np.median(par)
                sigma = np.std(par)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median, ),
                    r'$\sigma=%.2f$' % (sigma, )))
                ax.text(0.6, 0.85, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
                ax.set_xlabel(
                    r'$P(\frac{\theta_{obs}}{\theta_{0}})$', fontsize=18)
                for i in range(0, 3):
                    patches[i].set_facecolor('green')
                    patches[i].set_alpha(0.4)
                    patches[i].set_linewidth(1.5)
                    patches[i].set_edgecolor('green')
                for i in range(3, 11):
                    patches[i].set_facecolor('blue')
                    patches[i].set_alpha(0.4)
                    patches[i].set_linewidth(1.5)
                    patches[i].set_edgecolor('blue')
                for i in range(11, 40):
                    patches[i].set_facecolor('purple')
                    patches[i].set_alpha(0.4)
                    patches[i].set_linewidth(1.5)
                    patches[i].set_edgecolor('purple')
        if 'theta_obs' in choosen_pars:
            if T90_classification:
                s_GRB = np.asarray(s_GRB)
                l_GRB = np.asarray(l_GRB)
                w = 0.05
                n_s = math.ceil((s_GRB.max() - s_GRB.min())/w)
                n_l = math.ceil((l_GRB.max() - l_GRB.min())/w)
                weight_s = 1/(2*np.pi*np.sin(theta_obs_s))
                weight_l = 1/(2*np.pi*np.sin(theta_obs_l))
                plt.hist(s_GRB, bins=n_s, weights=weight_s,
                         histtype='step', linewidth=1.5, color='blue', alpha=0.7, edgecolor='blue')
                plt.hist(l_GRB, bins=n_l,
                         weights=weight_l, histtype='step', linewidth=1.5, color='blue', alpha=0.7, edgecolor='red')
                median_s = np.median(s_GRB)
                sigma_s = np.std(s_GRB)
                median_l = np.median(l_GRB)
                sigma_l = np.std(l_GRB)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_s, ),
                    r'$\sigma=%.2f$' % (sigma_s, )))
                textstr_l = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_l, ),
                    r'$\sigma=%.2f$' % (sigma_l, )))
                ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.6)
                ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.6)
            else:
                weight = 1/(2*np.pi*np.sin(theta_obs))
                plt.hist(par, bins=num_bins, color='blue',
                         weights=weight, alpha=0.5, edgecolor='blue', histtype='stepfilled', linewidth=1.5)
                plt.title('Probability to observe a certain theta_obs')
                median = np.median(par)
                sigma = np.std(par)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median, ),
                    r'$\sigma=%.2f$' % (sigma, )))
                ax.text(0.6, 0.85, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
                ax.set_xlabel(
                    r'$P(\theta_{obs})$', fontsize=18)
    if no_classification:
        '''
        if choosen_pars in Info['Log']:
            if T90_classification:
                plt.hist(pow(10, np.array(s_GRB)), bins=num_bins_s,
                         histtype='step', color='blue', alpha=0.7)
                plt.hist(pow(10, np.array(l_GRB)), bins=num_bins_l,
                         histtype='step', color='red', alpha=0.7)
                median_s = np.median(s_GRB)
                sigma_s = np.std(s_GRB)
                median_l = np.median(l_GRB)
                sigma_l = np.std(l_GRB)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_s, ),
                    r'$\sigma=%.2f$' % (sigma_s, )))
                textstr_l = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_l, ),
                    r'$\sigma=%.2f$' % (sigma_l, )))
                ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            else:
                plt.hist(pow(10, par), bins=num_bins, color='blue', alpha=0.4)
                median = np.median(par)
                sigma = np.std(par)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median, ),
                    r'$\sigma=%.2f$' % (sigma, )))
                ax.text(0.10, 0.85, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
        '''
        # else:
        if T90_classification:
            s_GRB = np.asarray(s_GRB)
            l_GRB = np.asarray(l_GRB)
            if choosen_pars == 'theta_obs' or choosen_pars == 'thetaobs_theta0' or choosen_pars == 'theta_0':
                w = 0.05
                n_s = math.ceil((s_GRB.max() - s_GRB.min())/w)
                n_l = math.ceil((l_GRB.max() - l_GRB.min())/w)
                plt.hist(s_GRB, bins=n_s,
                         histtype='step', linewidth=1.5, color='blue', alpha=0.7, edgecolor='blue')
                plt.hist(l_GRB, bins=n_l,
                         histtype='step', linewidth=1.5, color='red', alpha=0.7, edgecolor='red')
                median_s = np.median(s_GRB)
                sigma_s = np.std(s_GRB)
                median_l = np.median(l_GRB)
                sigma_l = np.std(l_GRB)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_s, ),
                    r'$\sigma=%.2f$' % (sigma_s, )))
                textstr_l = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_l, ),
                    r'$\sigma=%.2f$' % (sigma_l, )))
                ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                ax.set_xlabel(choosen_pars)
            else:
                w = 0.5
                n_s = math.ceil((s_GRB.max() - s_GRB.min())/w)
                n_l = math.ceil((l_GRB.max() - l_GRB.min())/w)
                plt.hist(s_GRB, bins=n_s,
                         histtype='step', linewidth=1.5, color='blue', alpha=0.7, edgecolor='blue')
                plt.hist(l_GRB, bins=n_l,
                         histtype='step', linewidth=1.5, color='red', alpha=0.7, edgecolor='red')
                median_s = np.median(s_GRB)
                sigma_s = np.std(s_GRB)
                median_l = np.median(l_GRB)
                sigma_l = np.std(l_GRB)
                textstr_s = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_s, ),
                    r'$\sigma=%.2f$' % (sigma_s, )))
                textstr_l = '\n'.join((
                    r'$\mathrm{median}=%.2f$' % (median_l, ),
                    r'$\sigma=%.2f$' % (sigma_l, )))
                # ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                # fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                # ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
                # fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                if 'E' in choosen_pars and 'Eta0' not in choosen_pars:
                    ax.set_xlabel(r'$\log(E_{50})$', fontsize=13)
                    ax.text(0.05, 0.98, textstr_s, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                    ax.text(0.05, 0.88, textstr_l, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                elif 'epsb' in choosen_pars:
                    ax.set_xlabel(r'$\log(\epsilon_{B})$', fontsize=13)
                    ax.text(0.60, 0.98, textstr_s, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                    ax.text(0.60, 0.88, textstr_l, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                elif 'epse' in choosen_pars:
                    ax.set_xlabel(r'$\log(\epsilon_{e})$', fontsize=13)
                    ax.text(0.60, 0.98, textstr_s, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                    ax.text(0.60, 0.88, textstr_l, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                elif 'n' in choosen_pars:
                    ax.set_xlabel(r'$\log(n_{0})$', fontsize=13)
                    ax.text(0.05, 0.98, textstr_s, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                    ax.text(0.05, 0.88, textstr_l, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                else:
                    ax.text(0.69, 0.98, textstr_s, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                    ax.text(0.69, 0.88, textstr_l, transform=ax.transAxes,
                            fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                    ax.set_xlabel(choosen_pars)
                # ax.set_xlabel(choosen_pars)

        else:
            plt.hist(par, bins=num_bins, color='blue',
                     alpha=0.5, edgecolor='blue', histtype='stepfilled', linewidth=1.5)
            median = np.median(par)
            sigma = np.std(par)
            textstr_s = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median, ),
                r'$\sigma=%.2f$' % (sigma, )))
            if 'E' in choosen_pars and 'Eta0' not in choosen_pars:
                ax.set_xlabel(r'$\log(E_{50})$', fontsize=13)
                ax.text(0.05, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
            elif 'epsb' in choosen_pars:
                ax.set_xlabel(r'$\log(\epsilon_{B})$', fontsize=13)
                ax.text(0.60, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
            elif 'epse' in choosen_pars:
                ax.set_xlabel(r'$\log(\epsilon_{e})$', fontsize=13)
                ax.text(0.60, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
            elif 'n' in choosen_pars:
                ax.set_xlabel(r'$\log(n_{0})$', fontsize=13)
                ax.text(0.05, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
            else:
                ax.text(0.69, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=15, verticalalignment='top', color='black')
                ax.set_xlabel(choosen_pars)

    if class_by_tSNE:
        plt.hist(class_small, bins=num_bins_class_small,
                 histtype='step', color='red', alpha=0.7)
        plt.hist(class_big, bins=num_bins_class_big,
                 histtype='step', color='blue', alpha=0.7)
        median_s = np.median(class_small)
        sigma_s = np.std(class_small)
        median_b = np.median(class_big)
        sigma_b = np.std(class_big)
        textstr_s = '\n'.join((
            r'$\mathrm{median}=%.2f$' % (median_s, ),
            r'$\sigma=%.2f$' % (sigma_s, )))
        textstr_b = '\n'.join((
            r'$\mathrm{median}=%.2f$' % (median_b, ),
            r'$\sigma=%.2f$' % (sigma_b, )))
        ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                fontsize=13, verticalalignment='top', color='red', alpha=0.7)
        ax.text(0.70, 0.89, textstr_b, transform=ax.transAxes,
                fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
    if class_by_gammaB:
        reasonable_short = np.asarray(reasonable_short)
        reasonable_long = np.asarray(reasonable_long)
        if choosen_pars == 'theta_obs' or choosen_pars == 'thetaobs_theta0' or choosen_pars == 'theta_0':
            w = 0.05
            n_s = math.ceil(
                (reasonable_short.max() - reasonable_short.min())/w)
            n_l = math.ceil((reasonable_long.max() - reasonable_long.min())/w)
            plt.hist(reasonable_short, bins=n_s,
                     histtype='step', linewidth=1.5, color='blue', alpha=0.7, edgecolor='blue')
            plt.hist(reasonable_long, bins=n_l,
                     histtype='step', linewidth=1.5, color='red', alpha=0.7, edgecolor='red')
            median_s = np.median(reasonable_short)
            sigma_s = np.std(reasonable_short)
            median_l = np.median(reasonable_long)
            sigma_l = np.std(reasonable_long)
            textstr_s = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_s, ),
                r'$\sigma=%.2f$' % (sigma_s, )))
            textstr_l = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_l, ),
                r'$\sigma=%.2f$' % (sigma_l, )))
            ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
            ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
                    fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            ax.set_xlabel(choosen_pars)
        else:
            w = 0.8
            n_s = math.ceil(
                (reasonable_short.max() - reasonable_short.min())/w)
            n_l = math.ceil((reasonable_long.max() - reasonable_long.min())/w)
            plt.hist(reasonable_short, bins=n_s,
                     histtype='step', linewidth=1.5, color='blue', alpha=0.7, edgecolor='blue')
            plt.hist(reasonable_long, bins=n_l,
                     histtype='step', linewidth=1.5, color='red', alpha=0.7, edgecolor='red')
            median_s = np.median(reasonable_short)
            sigma_s = np.std(reasonable_short)
            median_l = np.median(reasonable_long)
            sigma_l = np.std(reasonable_long)
            textstr_s = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_s, ),
                r'$\sigma=%.2f$' % (sigma_s, )))
            textstr_l = '\n'.join((
                r'$\mathrm{median}=%.2f$' % (median_l, ),
                r'$\sigma=%.2f$' % (sigma_l, )))
            # ax.text(0.70, 0.99, textstr_s, transform=ax.transAxes,
            # fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
            # ax.text(0.70, 0.89, textstr_l, transform=ax.transAxes,
            # fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            if 'E' in choosen_pars and 'Eta0' not in choosen_pars:
                ax.set_xlabel(r'$\log(E_{50})$', fontsize=13)
                ax.text(0.05, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.05, 0.88, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            elif 'epsb' in choosen_pars:
                ax.set_xlabel(r'$\log(\epsilon_{B})$', fontsize=13)
                ax.text(0.60, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.60, 0.88, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            elif 'epse' in choosen_pars:
                ax.set_xlabel(r'$\log(\epsilon_{e})$', fontsize=13)
                ax.text(0.60, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.60, 0.88, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            elif 'n' in choosen_pars:
                ax.set_xlabel(r'$\log(n_{0})$', fontsize=13)
                ax.text(0.05, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.05, 0.88, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
            else:
                ax.text(0.69, 0.98, textstr_s, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='blue', alpha=0.7)
                ax.text(0.69, 0.88, textstr_l, transform=ax.transAxes,
                        fontsize=13, verticalalignment='top', color='red', alpha=0.7)
                ax.set_xlabel(choosen_pars)
            # ax.set_xlabel(choosen_pars)
    return fig


def main():

    # config streamlit
    st.set_page_config(layout="wide")
    st.write("""
    # Fit parameters distributions
    """)

    choosen_pars = st.sidebar.selectbox(
        'Choose a parameter', ('redshift', 'Eta0', 'GammaB', 'theta_0', 'theta_obs', 'E', 'epsb', 'epse', 'n', 'p', 'thetaobs_theta0', 'GAMMA', 'chi2_red'))
    no_classification = st.sidebar.checkbox(
        'Display original distributions', True)
    T90_classification = st.sidebar.checkbox(
        'Long and short classification', True)
    prob_dist = st.sidebar.checkbox(
        'Display probability for theta_obs and theta_obs/theta_0', True)
    histo_classification = st.sidebar.checkbox(
        'Display distributions based on theta_obs/theta_0 classification', True)
    class_by_tSNE = st.sidebar.checkbox(
        'Display distributions with classification based on tSNE', True)
    class_by_gammaB = st.sidebar.checkbox(
        'Display distributions with classification based on gammaB >= 1/theta_obs criterion', True)

    # read data
    path_dirs = '/Users/alessandraberretta/'
    dirs_list = [path_dirs +
                 dir for dir in os.listdir(path_dirs) if dir.startswith('20')]
    summary_list = []
    DF1 = read_data_file(
        '/Users/alessandraberretta/JetFit/common_grb.csv')
    names_fit = DF1['common_grb'].values
    names_dist = []
    for elm in dirs_list:
        # if elm.endswith('results'):
        if 'results' in elm:
            for file in [elm + '/' + file for file in os.listdir(elm) if file.startswith('new_sum_') and file.endswith('.csv')]:
                summary_list.append(file)
                # names_dist.append(file.split('_')[6][0:-4])
                # print(file.split('_')[6][0:-4])
        # elif elm.endswith('simbad_res'):
        if 'simbad_res' in elm:
            for file in [elm + '/' + file for file in os.listdir(elm) if file.startswith('new_sum_') and file.endswith('.csv')]:
                summary_list.append(file)
                # names_dist.append(file.split('_')[7][0:-4])
                # print(file.split('_')[7][0:-4])

    for idx, elm in enumerate(summary_list):

        DF = read_data_file(elm)
        names_dist.append(elm.split('/')[4][12:-4])
        E.append(float(DF['Values'][3]))
        Eta0.append(float(DF['Values'][0]))
        GammaB.append(float(DF['Values'][1]))
        epsB.append(float(DF['Values'][4]))
        epse.append(float(DF['Values'][5]))
        n.append(float(DF['Values'][6]))
        p.append(float(DF['Values'][7]))
        theta_obs.append(float(DF['Values'][2]))
        redshift.append(float(DF['Values'][8]))
        chi2_red.append(float(DF['Values'][12]))
        Class_T90.append(DF['Values'].iloc[-1])

    theta_0 = []
    thetaobs_theta0 = []
    thetaobs_theta0_rad = []
    GAMMA = []

    theta_obs_deg = np.rad2deg(theta_obs)

    for gamma in GammaB:
        theta_0.append(1/gamma)

    theta_0_deg = np.rad2deg(theta_0)

    for idx, par in enumerate(theta_obs_deg):
        thetaobs_theta0.append(np.float64(
            par)/(np.float64(theta_0_deg[idx]))/2)

    for idx, par in enumerate(theta_obs):
        thetaobs_theta0_rad.append(np.float64(
            par)/(np.float64(theta_0[idx]))/2)

    for idx_eta, eta in enumerate(Eta0):
        GAMMA.append(2*np.float64(eta)*np.float64(GammaB[idx]))

    Pars = {
        'E': np.array(E),
        'n': np.array(n),
        'Eta0': np.array(Eta0),
        'GammaB': np.array(GammaB),
        'theta_0': np.array(theta_0),
        'theta_obs': np.array(theta_obs),
        'epse': np.array(epse),
        'epsb': np.array(epsB),
        'p': np.array(p),
        'thetaobs_theta0': np.array(thetaobs_theta0),
        'GAMMA': np.array(GAMMA),
        'chi2_red': np.array(chi2_red),
        'Class_T90': np.array(Class_T90),
        'redshift': np.array(redshift),
        'names': np.array(names_dist)
    }

    print('tot', len(names_dist))

    df = pd.DataFrame()
    df['names'] = names_dist
    df.to_csv('prova_2.csv', sep='\t')

    inf_in_chi2red = np.isinf(Pars['chi2_red'])

    for idx, elm in enumerate(inf_in_chi2red):
        if elm == True:
            Pars['chi2_red'][idx] = -2

    conta = []
    for idx, elm in enumerate(thetaobs_theta0):
        if elm > 1:
            conta.append(redshift[idx])
    # print(conta)

    l_GRB = []
    s_GRB = []
    theta_obs_s = []
    theta_obs_l = []
    class_1_temp = []
    class_2_temp = []
    class_small = []
    class_big = []
    reasonable_short = []
    reasonable_long = []

    for idx, elm in enumerate(Pars['Class_T90']):
        if elm == 'S':
            s_GRB.append(Pars[choosen_pars][idx])
            if (Pars['theta_obs'][idx]/(Pars['theta_0'][idx])/2) > 1:
                print('s', Pars['names'][idx])
            elif Pars['GammaB'][idx] <= 1/(Pars['theta_obs'][idx]*2):
                reasonable_short.append(Pars[choosen_pars][idx])
        else:
            l_GRB.append(Pars[choosen_pars][idx])
            if Pars['GammaB'][idx] <= 1/(Pars['theta_obs'][idx]*2):
                reasonable_long.append(Pars[choosen_pars][idx])

    print('short', len(reasonable_long))
    print('long', len(reasonable_short))

    for idx, elm in enumerate(Pars['Class_T90']):
        if elm == 'S':
            theta_obs_s.append(Pars['theta_obs'][idx])
        else:
            theta_obs_l.append(Pars['theta_obs'][idx])

    for idx, elm in enumerate(Pars['thetaobs_theta0']):
        if elm > 0.21:
            class_1_temp.append(Pars[choosen_pars][idx])
        else:
            class_2_temp.append(Pars[choosen_pars][idx])

    for idx, elm in enumerate(Pars['names']):
        if elm == '131105A' or elm == '050318' or elm == '091127' or elm == '150821A' or elm == '120422A' or elm == '140903A' or elm == '130603B' or elm == '081203A' or elm == '090510' or elm == '060604' or elm == '101219A' or elm == '120909A':
            class_small.append(Pars[choosen_pars][idx])
        else:
            class_big.append(Pars[choosen_pars][idx])

    Latex = {
        'E': r'$log_{10}E_{j,50}$',
        'n': r'$log_{10}n_{0,0}$',
        'Eta0': r'$\eta_0$',
        'GammaB': r'$\gamma_B$',
        'theta_0': r'$\theta_{0}$',
        'theta_obs': r'$\theta_{obs}$',
        'epse': r'$log_{10}\epsilon_e$',
        'epsb': r'$log_{10}\epsilon_B$',
        'p': r'$p$',
        'thetaobs_theta0': r'$frac{\theta_{obs}}{\theta_{0}}$',
        'GAMMA': r'$\Gamma$',
        'chi2_red': r'$\chi^{2}_{red}$',
        'Class_T90': r'$T_{90}$',
        'redshift': r'$z$'
    }
    # num_bins_s=n_s, num_bins_l=n_l,
    Label = []
    for x in Info['Fit']:
        Label.append(Latex[x])

    hist = hist_pars(class_by_tSNE, class_by_gammaB, no_classification, class_1_temp, class_2_temp, prob_dist, histo_classification,
                     class_small, class_big, Pars[choosen_pars], Info, choosen_pars, theta_obs, theta_obs_s,
                     theta_obs_l, s_GRB, l_GRB, reasonable_long, reasonable_short, T90_classification, num_bins=40, num_bins_class_small=30, num_bins_class_big=30)

    st.pyplot(hist)


if __name__ == "__main__":
    main()
