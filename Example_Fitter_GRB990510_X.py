import corner
import numpy as np
import emcee as em
import matplotlib.pyplot as plt
import pandas as pd
import sys

from JetFit import FitterClass

# Parameters
Table = "./Table.h5"

Info = {
    # Fitting parameters (Parameter names see P dictionary below)
    'Fit': np.array(['Eta0', 'GammaB', 'theta_obs']),
    # Set parameters in log scale
    'Log': np.array(['E', 'n', 'epse', 'epsb']),
    'LogType': 'Log10',                              # Log scale type: Log10 or Log
    # Prior for observation angle: Sine or Uniform
    'ThetaObsPrior': 'Sine',
    # Flux type: Spectral or Integrated
    'FluxType': 'Spectral'
}

# Bounds for parameters. All in linear scale.
FitBound = {
    'E': np.array([1e-6, 1e3]),
    'n': np.array([1e-6, 1e3]),
    'Eta0': np.array([2., 10.]),
    'GammaB': np.array([1., 12.]),
    'theta_obs': np.array([0., 1]),
    'epse': np.array([1e-6, 1.]),
    'epsb': np.array([1e-6, 1.]),
    'p': np.array([2., 4.])
}

# P:
# For non-fiting parameters, P set default values.
# For fitting paramters, P:
#  1. If Explore == True: Fitting parameters are randomly distributed in whole parameter space.
#  2. If Explore != True: Fitting parameters are randomly distributed around maximum posterior region, indicated by values in P.

Explore = True

# E_loop = np.linspace(0.15, 5.0, 10)
# Eta0_loop = np.linspace(2., 10., 9)
# GammaB_loop = np.linspace(1., 12., 12)

# P =

# for E in E_loop:
# for Eta0 in Eta0_loop:
# for GammaB in GammaB_loop:

P = {
    'E': 1.15869069395227384,
    'Eta0': 9.973477192135503,
    'GammaB': 5.000923300022666,
    'dL': 3.761,
    'epsb': 0.8,
    'epse': 0.8,
    'n': 0.08,
    'p': 2.2,
    'theta_obs': 0.005,
    'xiN': 1.0,
    'z': 1.619
}

# P.append(P_temp)

# for p in P:

# parameters for fitter.
# Path to observation data.
GRB = './GRB990510_X.csv'

# for demostaration
SamplerType = "ParallelTempered"
NTemps = 2
NWalkers = 6
Threads = 8

BurnLength = 10
RunLength = 10

# For GW170817 and bellow parameters, it takes ~24 hours to finish.
# For quick run, values of the parameters can be modified accordingly.
# SamplerType = "ParallelTempered"
# NTemps = 10
# NWalkers = 100
# Threads = 8

# BurnLength = 10000
# RunLength = 10000

# Fitter
# Initialize Fitter
Fitter = FitterClass(Table, Info, FitBound, P, Explore=Explore)
# LoadData
DF = pd.read_csv(GRB)
Times, Fluxes, FluxErrs, Freqs = DF['Times'].values, DF['Fluxes'].values * \
    1e3, DF['FluxErrs'].values*1e3, DF['Freqs'].values
Fitter.LoadData(Times, Fluxes, FluxErrs, Freqs)
# Initialize sampler
Fitter.GetSampler(SamplerType, NTemps, NWalkers, Threads)

# Fitting
# Burning in
BurnInResult = Fitter.BurnIn(BurnLength=BurnLength)
# Fitting and store chain results to Result
Result = Fitter.RunSampler(RunLength=RunLength, Output=None)

# Analysis: Below only depends on Result
# Plot Light Curves


def Log2Linear(Log, Info):
    Linear = []
    for i, key in enumerate(Info['Fit']):
        if key in Info['Log']:
            if Info['LogType'] == 'Log10':
                Linear.append(np.power(10., Log[i]))
            else:
                Linear.append(np.exp(Log[i]))
        else:
            Linear.append(Log[i])
    return np.array(Linear)


def PltDF(ax, DF, ColorList=['orange', 'red', 'g', 'b', 'black', 'pink'], ScaleFactor=[1., 1., 1., 1., 1., 1.], Legend=True, XAxisDay=False):
    Freqs = DF['Freqs'].unique()

    for Freq, Color, Scale in zip(Freqs, ColorList, ScaleFactor):
        SubDF = DF[DF['Freqs'] == Freq]
        if XAxisDay:
            Times = SubDF['Times']/24./3600
        else:
            Times = SubDF['Times']
        Fluxes = SubDF['Fluxes']
        FluxErrs = SubDF['FluxErrs']
        if max(ScaleFactor) > 1.:
            label = '%.1e x %d' % (Freq, Scale)
        else:
            label = '%.1e' % Freq
        ax.errorbar(Times, Fluxes*Scale, yerr=FluxErrs *
                    Scale, color=Color, fmt='.', label=label)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel('Flux density (mJy)')
    # ax.set_ylim(1e-7, 1e-3)
    if XAxisDay:
        ax.set_xlabel('Time (day)')
    else:
        ax.set_xlabel('Time (s)')
    if Legend:
        ax.legend(loc=0)


# Find the best fitting parameters
TheChain = Result['Chain']
LnProbability = Result['LnProbability']
FitDim = len(Info['Fit'])

BestWalker = np.unravel_index(
    np.nanargmax(LnProbability), LnProbability.shape)
BestParameter = TheChain[BestWalker]
BestLnProbability = LnProbability[BestWalker]
BestLinearParameter = Log2Linear(BestParameter, Info)

BestP = P.copy()
for i, key in enumerate(Info['Fit']):
    BestP[key] = BestLinearParameter[i]

# Plot best fitting light curves
fig, ax = plt.subplots(figsize=(8, 8))
# ColorList = ['orange', 'red', 'g', 'b']
# ScaleFactor = [6., 1., 100., 800.]
ColorList = ['orange']
ScaleFactor = [1.]

PltDF(ax, DF, ColorList=ColorList,
      ScaleFactor=ScaleFactor, Legend=True, XAxisDay=False)

NPoints = 200
Left = 1.
Right = 2.
for i, Freq in enumerate(DF['Freqs'].unique()):
    idx = DF['Freqs'] == Freq
    NewTimes = np.linspace(DF['Times'].min()*Left,
                           DF['Times'].max()*Right, NPoints)
    NewFreqs = np.ones(len(NewTimes))*Freq

    # Generate Fluxes
    FluxesModel = np.asarray(
        Fitter.FluxGenerator.GetSpectral(NewTimes, NewFreqs, BestP))
    print(FluxesModel)

    plt.loglog(NewTimes, FluxesModel *
               ScaleFactor[i], '--', color=ColorList[i], linewidth=1.5)

plt.savefig('light_curves_GRB990510_X.png')
# plt.savefig('/Users/alessandraberretta/Desktop/JetFit/990510/light_curves_GRB990510_{}_{}_{}.png'.format(
# p['E'], p['Eta0'], p['GammaB']))

# Plot Distribution
# Get nice latex label
Latex = {
    'E': r'$log_{10}E_{j,50}$',
    'n': r'$log_{10}n_{0,0}$',
    'Eta0': r'$\eta_0$',
    'GammaB': r'$\gamma_B$',
    'theta_obs': r'$\theta_{obs}$',
    'epse': r'$log_{10}\epsilon_e$',
    'epsb': r'$log_{10}\epsilon_B$',
    'p': r'$p$'
}

Label = []
for x in Info['Fit']:
    Label.append(Latex[x])

# plot contour with ChainConsumer
Chain = Result['Chain'].reshape((-1, FitDim))
fig = corner.corner(Chain, labels=Label, label_size=20, bins=40, plot_datapoints=False,
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, color='darkblue',
                    label_kwargs={'fontsize': 18},
                    title_kwargs={"fontsize": 18})

fig.savefig('contour_GRB990510_X.png')
# fig.savefig('/Users/alessandraberretta/Desktop/JetFit/990510/contour_GRB990510_{}_{}_{}.png'.format(
# p['E'], p['Eta0'], p['GammaB']))
