import sys
import numpy as np
from scipy.stats import linregress
import statistics
import matplotlib.pyplot as plt
from plot_analysis import lc_plot


# function to compute the standard deviation for bunch of points
def stdev_groups(GRB, times_mean, fluxes_mean, fluxerrs_mean, group_points):

    print(len(times_mean))
    '''
    if len(times_mean) % group_points != 0:
        st.error("change the number of points to group with")
        sys.exit()
    '''

    slopes = []
    intercepts = []
    # times_rebin = []
    first_point_t = []
    last_point_t = []
    first_point_f = []
    last_point_f = []
    first_flux_errs = []
    last_flux_errs = []

    for idx in range(len(times_mean)-1):
        slope, intercept, _, _, _ = linregress([np.log10(times_mean[idx]), np.log10(
            times_mean[idx+1])], [np.log10(fluxes_mean[idx]), np.log10(fluxes_mean[idx+1])])
        slopes.append(slope)
        intercepts.append(intercept)
        # times_rebin.append(times_mean[idx])
        first_point_t.append(times_mean[idx])
        last_point_t.append(times_mean[idx+1])
        first_point_f.append(fluxes_mean[idx])
        last_point_f.append(fluxes_mean[idx+1])
        first_flux_errs.append(fluxerrs_mean[idx])
        last_flux_errs.append(fluxerrs_mean[idx+1])

    grouped_slopes = [slopes[n:n+group_points]
                      for n in range(0, len(slopes), group_points)]
    # grouped_times_rebin = [times_rebin[n:n+group_points]
    # for n in range(0, len(times_rebin), group_points)]
    grouped_first_point_t = [first_point_t[n:n+group_points]
                             for n in range(0, len(first_point_t), group_points)]
    grouped_last_point_t = [last_point_t[n:n+group_points]
                            for n in range(0, len(last_point_t), group_points)]
    grouped_first_point_f = [first_point_f[n:n+group_points]
                             for n in range(0, len(first_point_f), group_points)]
    grouped_last_point_f = [last_point_f[n:n+group_points]
                            for n in range(0, len(last_point_f), group_points)]
    grouped_first_fluxerrs = [first_flux_errs[n:n+group_points]
                              for n in range(0, len(first_flux_errs), group_points)]
    grouped_last_fluxerrs = [last_flux_errs[n:n+group_points]
                             for n in range(0, len(last_flux_errs), group_points)]

    stdev_grouped = []
    for _list in grouped_slopes:
        stdev_grouped.append(statistics.stdev(_list))

    mean_stdev_grouped = sum(stdev_grouped) / len(stdev_grouped)
    dropped_idx = []
    '''
    good_slopes = np.empty(0)

    # print(grouped_times_rebin)

    for idx, elm in enumerate(grouped_slopes):
        if idx != 1:
            good_slopes = np.append(good_slopes, np.asarray(elm))

    std = np.std(good_slopes)
    # print(2*std)

    slopes_arr = np.array(slopes)
    mask_slopes = np.abs(slopes_arr) < (2*std)
    good_slopes_real = slopes_arr[mask_slopes]

    good_times_rebin = []
    good_first_point_f = np.empty(0)
    good_first_fluxerrs = np.empty(0)

    for _list in grouped_times_rebin:
        temp = np.asarray(_list)

    t_rebin_arr = np.asarray(grouped_times_rebin)
    t_rebin_flatten = t_rebin_arr.flatten()
    print(t_rebin_flatten)

    for idx, elm in enumerate(good_slopes_real):
        print('ciao')
        good_times_rebin.append(t_rebin_flatten[idx])
    # print(good_times_rebin)
    '''
    for elm in stdev_grouped:
        if elm > mean_stdev_grouped + np.std(stdev_grouped):
            dropped_idx.append(stdev_grouped.index(elm))

    grouped_slopes = np.delete(grouped_slopes, dropped_idx)

    grouped_slopes = np.hstack(grouped_slopes)
    # grouped_times_rebin = np.hstack(grouped_times_rebin)
    grouped_first_point_t = np.hstack(grouped_first_point_t)
    grouped_last_point_t = np.hstack(grouped_last_point_t)
    grouped_first_point_f = np.hstack(grouped_first_point_f)
    grouped_last_point_f = np.hstack(grouped_last_point_f)
    grouped_first_fluxerrs = np.hstack(grouped_first_fluxerrs)
    grouped_last_fluxerrs = np.hstack(grouped_last_fluxerrs)
    std_grouped_slopes_removed = np.std(grouped_slopes)
    print(std_grouped_slopes_removed)

    # slopes_red_bystd = []
    # dropped_idx_bystd = []
    Times_red_bystd = []
    Fluxes_red_bystd = []
    FluxErrs_red_bystd = []
    print(times_mean)
    print(grouped_first_point_t)
    print(grouped_last_point_t)
    print(slopes)
    for idx, elm in enumerate(slopes):
        if abs(elm) > std_grouped_slopes_removed:
            if elm > 0:
                Times_red_bystd.append(grouped_first_point_t[idx])
                Fluxes_red_bystd.append(grouped_first_point_f[idx])
                FluxErrs_red_bystd.append(grouped_first_fluxerrs[idx])
            else:
                Times_red_bystd.append(grouped_last_point_t[idx])
                Fluxes_red_bystd.append(grouped_last_point_f[idx])
                FluxErrs_red_bystd.append(grouped_last_fluxerrs[idx])
            # slopes_red_bystd.append(elm)
            # dropped_idx_bystd.append(idx)

    # Times_red_bystd = np.delete(grouped_first_point_t, dropped_idx_bystd)
    # Fluxes_red_bystd = np.delete(grouped_first_point_f, dropped_idx_bystd)
    # FluxErrs_red_bystd = np.delete(grouped_first_fluxerrs, dropped_idx_bystd)

    '''
    slopes_rem = []
    times_rebin_rem = []
    first_point_t_rem = []
    last_point_t_rem = []
    first_point_f_rem = []
    first_fluxerrs_rem = []
    last_fluxerrs_rem = []
    for idx, _ in enumerate(grouped_slopes):
        slopes_rem += grouped_slopes[idx]
        times_rebin_rem += grouped_times_rebin[idx]
        first_point_t_rem += grouped_first_point_t[idx]
        last_point_t_rem += grouped_last_point_t[idx]
        first_point_f_rem += grouped_first_point_f[idx]
        first_fluxerrs_rem += grouped_first_fluxerrs[idx]
        last_fluxerrs_rem += grouped_last_fluxerrs[idx]
    '''

    lc_grouped = lc_plot(GRB, Times_red_bystd, Fluxes_red_bystd, FluxErrs_red_bystd,
                         'blue', original=False, rebin=False, flare=False, group=True)

    return stdev_grouped, mean_stdev_grouped, std_grouped_slopes_removed, lc_grouped
