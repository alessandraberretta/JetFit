import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from decimal import Decimal


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
