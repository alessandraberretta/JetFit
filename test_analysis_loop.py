from random import randint
import corner
import numpy as np
import emcee as em
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math
from numpy import random
from decimal import Decimal
from JetFit import FitterClass
import itertools

# Parameters
Table = "./Table.h5"

Info = {
    # Fitting parameters (Parameter names see P dictionary below)
    'Fit': np.array(['Eta0', 'GammaB', 'theta_obs', 'E', 'epsb', 'epse', 'n', 'p']),
    # Set parameters in log scale
    'Log': np.array(['E', 'n', 'epse', 'epsb']),
    'LogType': 'Log10',                              # Log scale type: Log10 or Log
    # Prior for observation angle: Sine or Uniform
    'ThetaObsPrior': 'Uniform',
    # Flux type: Spectral or Integrated
    'FluxType': 'Integrated'
}

# Bounds for parameters. All in linear scale.
FitBound = {
    'E': np.array([1e-8, 1e3]),
    'n': np.array([1e-5, 1e4]),
    'Eta0': np.array([2., 10.]),
    'GammaB': np.array([1., 12.]),
    'theta_obs': np.array([0.045, 1.]),
    'epse': np.array([1e-6, 50.]),
    'epsb': np.array([1e-6, 50.]),
    'p': np.array([2., 5.])
}


def percentage(num, val_perc):
    return float(num) * float(val_perc)/100


def random_generator(min, max, num, seed=None):
    points = []
    np.random.seed(seed)
    while len(points) != num:
        val = np.random.randint(min, max)
        if val not in points:
            points.append(val)

    return np.asarray(points)


# P:
# For non-fiting parameters, P set default values.
# For fitting paramters, P:
#  1. If Explore == True: Fitting parameters are randomly distributed in whole parameter space.
#  2. If Explore != True: Fitting parameters are randomly distributed around maximum posterior region, indicated by values in P.
# p_loop = np.random.normal(2.6333493591554804, 0.52666987183, 500)
# theta_obs_loop = np.random.normal(0.04769798916899842, 0.00953959783, 500)

Explore = False

E_loop = np.linspace(0.09521441637, 0.22216697153, 100)
Eta0_loop = np.linspace(5.81878175371, 8.72817263057, 100)
GammaB_loop = np.linspace(5.60073864002, 8.40110796002, 100)
epsb_loop = np.linspace(0.00106589652, 0.00159884478, 100)
epse_loop = np.linspace(0.02443670305, 0.05701897379, 100)
n_loop = np.linspace(0.07896976823, 0.11845465234, 100)
p_loop = np.linspace(2.18667948733, 3.28001923099, 100)
theta_obs_loop = np.linspace(0.014158391336, 0.02123758699, 100)


P_random = []

P1 = {
    'E': 0.15869069395227384,
    'Eta0': 7.273477192135503,
    'GammaB': 7.000923300022666,
    'dL': 0.00264908,
    'epsb': 0.0013323706571267526,
    'epse': 0.04072783842837688,
    'n': 0.09871221028954489,
    'p': 2.7333493591554804,
    'theta_obs': 0.01769798916899842,
    'xiN': 1.0,
    'z': 0.002
}

P2 = {
    'E': 0.015869069395227384,
    'Eta0': 5.973477192135503,
    'GammaB': 9.000923300022666,
    'dL': 0.00264908,
    'epsb': 0.0013323706571267526,
    'epse': 0.05072783842837688,
    'n': 1.9871221028954489,
    'p': 2.4333493591554804,
    'theta_obs': 0.1,
    'xiN': 1.0,
    'z': 0.002
}


for i in range(100):

    tmp_P = {
        'E': np.random.choice(E_loop),
        'Eta0': np.random.choice(Eta0_loop),
        'GammaB': np.random.choice(GammaB_loop),
        'dL': 0.00264908,
        'epsb': np.random.choice(epsb_loop),
        'epse': np.random.choice(epse_loop),
        'n': np.random.choice(n_loop),
        'p': np.random.choice(p_loop),
        'theta_obs': np.random.choice(theta_obs_loop),
        'xiN': 1.0,
        'z': 0.002
    }

    P_random.append(tmp_P)

fig, ax = plt.subplots(figsize=(8, 8))
ColorList = ['black']
ScaleFactor = [1.]
start_point = 5.12e3
end_point = 6.78e6

NewTimes = np.arange(start_point, end_point,
                     step=(end_point-start_point)/10000)

Fitter = FitterClass(Table, Info, FitBound, P1, Explore=Explore)

for idx, P in enumerate(P_random):

    data_P = []
    for idx_par, par in enumerate(Info['Fit']):
        data_P.append([par, P[par]])
    df_P = pd.DataFrame(data_P, columns=['Parameters', 'Values'])
    df_P.to_csv('/Users/alessandraberretta/JetFit/data_100_all_P1/summary_P1/summary_P{}_all_P1.csv'.format(idx),
                index=False,  sep='\t')

    points = random_generator(10, 100, 1, seed=np.random.randint(1, 100000))

    for val in points:

        NewTimes_red = np.sort(
            np.random.choice(NewTimes, val))
        Freqs = np.ones(len(NewTimes_red))*2e18
        NewFreqs = np.repeat([[7.254E+16, 2.418E+18]],
                             len(NewTimes_red), axis=0)
        partial = Fitter.FluxGenerator.GetIntegratedFlux(
            NewTimes_red, NewFreqs, P)
        FluxesModel = np.asarray(partial)
        if np.isnan(FluxesModel).any():
            continue
        print(FluxesModel)
        FluxesModel_diff = []
        for x in FluxesModel:
            perc = percentage(x, 25)
            FluxesModel_diff.append(
                x+random.uniform(-perc, perc))
        FluxesModel_diff_err = np.asarray(
            FluxesModel_diff)
        plt.loglog(NewTimes_red, FluxesModel*ScaleFactor,
                   '--', color='black', linewidth=1.5)
        plt.scatter(NewTimes_red, FluxesModel, c='red')
        ax = plt.gca()
        ax.errorbar(NewTimes_red, FluxesModel_diff,
                    yerr=0.25*FluxesModel_diff_err, fmt='.', color='black')
        # plt.errorbar(NewTimes, FluxesModel_diff, yerr=FluxesModel_diff_err, fmt='o')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel(
            'Flux (0.3 - 10 keV) (erg/cm^2/s)')
        ax.set_xlabel('Time (s)')
        plt.savefig(
            '/Users/alessandraberretta/JetFit/data_100_all_P1/lc_P1/lc_P{}_{}_all_P1.png'.format(idx, val))
        plt.clf()
        d = {'Times': NewTimes_red, 'Freqs': Freqs,
             'Fluxes': FluxesModel_diff, 'FluxErrs': 0.25*FluxesModel_diff_err}
        df = pd.DataFrame(
            data=d, columns=['Times', 'Freqs', 'Fluxes', 'FluxErrs'])
        df["Times"] = ['%.6E' %
                       Decimal(x) for x in df['Times']]
        df["Freqs"] = ['%.6E' %
                       Decimal(x) for x in df['Freqs']]
        df["Fluxes"] = ['%.6E' %
                        Decimal(y) for y in df['Fluxes']]
        df["FluxErrs"] = ['%.6E' %
                          Decimal(z) for z in df['FluxErrs']]
        d2 = {'Times': NewTimes_red, 'FluxesModel': FluxesModel}
        df2 = pd.DataFrame(
            data=d2, columns=['Times', 'FluxesModel'])
        df2["Times"] = ['%.6E' %
                        Decimal(x) for x in df2['Times']]
        df2["FluxesModel"] = ['%.6E' %
                              Decimal(x) for x in df2['FluxesModel']]
        df.to_csv(
            '/Users/alessandraberretta/JetFit/data_100_all_P1/data_100_all_P1/data_P{}_{}_all_P1.csv'.format(idx, val), index=False)
        df2.to_csv(
            '/Users/alessandraberretta/JetFit/data_100_all_P1/data_100_all_P1/tf_P{}_{}_all_P1.csv'.format(idx, val), index=False)
