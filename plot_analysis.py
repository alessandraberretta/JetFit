import matplotlib.pyplot as plt
import numpy as np


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
def lc_plot(GRB, times, fluxes, fluxes_err, color, original=True, rebin=True, flare=True, group=True):
    fig, ax = plt.subplots(clear=True)
    # plt.clf()
    ax.errorbar(times, fluxes, yerr=fluxes_err, fmt='.', color=color)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('times (s)')
    ax.set_ylabel('flux (erg/cm^2/s)')
    if original:
        plt.title('Original lightcurve')
        # plt.savefig("lc_" + GRB + ".pdf")
    if rebin:
        plt.title('Re-binned lightcurve')
        # plt.savefig("lc_" + GRB + "_rebin.pdf")
    if flare:
        plt.title('Removed flares lightcurve')
        # plt.savefig("lc_" + GRB + "_nf.pdf")
    if group:
        # plt.xlim(1e2, 2*np.max(times))
        plt.title('Flare removal after slopes grouping')
    return fig
