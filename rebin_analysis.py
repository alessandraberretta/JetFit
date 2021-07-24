import numpy as np
import matplotlib.pyplot as plt
from plot_analysis import lc_plot


# function for weighted mean computation
def weighted_error_mean(var, err_var):
    return np.average(var, axis=0, weights=[1/i for i in err_var])


# function for data rebinning
def rebin_data(GRB, DF, t_0, t_rebin, partial_rebin):

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

    if t_rebin == 1:
        times_mean.append(Times[-1])
        fluxes_mean.append(Fluxes[-1])
        fluxerrs_mean.append(FluxErrs[-1])

    # get the rebinned lightcurve
    lc_rebin = lc_plot(GRB, times_mean, fluxes_mean,
                       fluxerrs_mean, 'red', original=False, rebin=True, flare=False, group=False)

    return lc_rebin, times_mean, fluxes_mean, fluxerrs_mean
