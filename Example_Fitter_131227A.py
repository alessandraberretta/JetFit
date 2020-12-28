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
    'Fit': np.array(['Eta0', 'GammaB', 'theta_obs', 'E', 'n', 'epse', 'epsb', 'p']),
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
    'epsb': np.array([1e-10, 1.]),
    'p': np.array([2., 4.])
}


# P:
# For non-fiting parameters, P set default values.
# For fitting paramters, P:
#  1. If Explore == True: Fitting parameters are randomly distributed in whole parameter space.
#  2. If Explore != True: Fitting parameters are randomly distributed around maximum posterior region, indicated by values in P.

Explore = True

P = {
    'E': 1.9869069395227384,
    'Eta0': 9.973477192135503,
    'GammaB': 11.000923300022666,
    'dL': 0.396680808,
    'epsb': 0.00065,
    'epse': 0.0065,
    'n': 1000,
    'p': 2.22,
    'theta_obs': 0.088,
    'xiN': 1.0,
    'z': 5.3
}

# parameters for fitter.
# Path to observation data.
GRB = './GRB131227A_new.csv'

# for demostaration
SamplerType = "ParallelTempered"
NTemps = 84
NWalkers = 284
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
Times, Fluxes, FluxErrs, Freqs = DF['Times'].values, DF['Fluxes'].values, DF['FluxErrs'].values, DF['Freqs'].values
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


def PltDF(ax, DF, ColorList=['orange', 'red', 'g', 'b'], ScaleFactor=[1., 1., 1., 1.], Legend=True, XAxisDay=False):
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

BestWalker = np.unravel_index(np.nanargmax(LnProbability), LnProbability.shape)
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
ColorList = ['purple']
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
    # print(FluxesModel)

    plt.loglog(NewTimes, FluxesModel *
               ScaleFactor[i], '--', color=ColorList[i], linewidth=1.5)

plt.savefig('light_curve_GRB131227A_fit.png')

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

fig.savefig("contour_GRB131227A_fit.png")
