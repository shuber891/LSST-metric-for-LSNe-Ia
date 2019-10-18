#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:30:24 2019

@author: shuber
"""

import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt

# To run this code you have to download the files from https://me.lsst.eu/gris/
# and chance the data_path and plot_path files, further you have to specifiy
# the area in the d_area_observing_strategies dic if area is different from
# the 18000 deg^2 nominal WFD area

# This is where your cadence strategy file is located
data_path = "/afs/mpa/data/shuber/data/LSST_MAF_data/"
plot_path = "/afs/mpa/home/shuber/research/MicrolensSN/plot/"

# list of observing stratigies to investigate where *.npy file
# has been stored in your data path
observing_strategies = ["alt_sched", "alt_sched_rolling", "kraken_2026",
                        "kraken_2035", "colossus_2667", "pontus_2002",
                        "pontus_2489", "mothra_2049", "nexus_2097",
                        "kraken_2044", "kraken_2042", "rolling_10yrs_opsim",
                        "rolling_mix_10yrs_opsim"]
observing_strategies = ["alt_sched", "alt_sched_rolling", "kraken_2026",
                        "colossus_2667", "rolling_10yrs_opsim"]

# This dictionary contains the WFD survey area for cadence strategies.
# If strategy is not in that dic the nominal value 18000 deg^2 will be assumed
d_area_observing_strategies = {"kraken_2044": 24700, "pontus_2002": 24700,
                               "nexus_2097": 24700, "mothra_2049": 24700}


def f_get_median(array):
    l_median = np.ndarray.tolist(array)
    l_median.sort()
    median = l_median[int(len(l_median)*0.5)] * array.unit
    return median


def f_fit_function(gap_median_filter):
    # has been constructed by fitting 13 observing strategies
    fit = 2.14 * np.exp(0.36 * gap_median_filter.value)
    return fit


d_good_delay_from_fit = {}
for observing_strategy in observing_strategies:
    if observing_strategy in d_area_observing_strategies.keys():
        survey_area = d_area_observing_strategies[observing_strategy] * u.deg**2
    else:
        survey_area = 18000 * u.deg**2
    array = np.load("%s%s_SLMetric_new.npy" % (data_path, observing_strategy))

    gap_median_filter = f_get_median(array["gap_median"] * u.day)
    cumulative_season_length = f_get_median((array["season_length"] * u.day).to(u.year))

    """
    #Further parameters stored in *.npy but not used
    area_pixel = array["area"] * u.deg**2
    healpixId = array["healpixId"]
    pixRa = array["pixRa"] * u.deg
    pixDec = array["pixDec"] *u.deg
    gap_max_filter = array["gap_max"] * u.day
    gap_median_u =  array["gap_median_u"] * u.day
    gap_median_g =  array["gap_median_g"] * u.day
    gap_median_r =  array["gap_median_r"] * u.day
    gap_median_i =  array["gap_median_i"] * u.day
    gap_median_z =  array["gap_median_z"] * u.day
    gap_median_y =  array["gap_median_y"] * u.day
    """

    N_OM10 = 45.7
    survey_area_OM10 = 20000 * u.deg**2
    cumulative_season_length_OM10 = 2.5 * u.year

    # rescale numbers of OM10 with cumulative_season_length and survey_area
    # of a given cadence strategy
    total_number_of_LSNeIa_0M10 = N_OM10 * survey_area/survey_area_OM10 * cumulative_season_length/cumulative_season_length_OM10

    # this estimates the number of well LSNeIa which have a well measured time
    # delay, i.e. that at least one image as accuracy < 1 percent
    # and precision < 5 percent
    number_of_good_delay_LSNeIa = total_number_of_LSNeIa_0M10 / f_fit_function(gap_median_filter)

    d_good_delay_from_fit[observing_strategy] = number_of_good_delay_LSNeIa
    # print observing_strategy, survey_area, cumulative_season_length, gap_median_filter, number_of_good_delay_LSNeIa

# This part creates a bar plot
observing_strategies_ord = []
y_follow_up = []
x_axis = []
x_ticks = []
x_counter = 0
for observing_strategy, total_number in sorted(d_good_delay_from_fit.iteritems(),reverse=True,key=lambda (k,v): (v,k)):  
    observing_strategies_ord.append(observing_strategy)

    x_counter += 2
    x_axis.append(x_counter)
    x_ticks.append(observing_strategy)
    y_follow_up.append(total_number)

x_axis = np.asarray(x_axis)
ticks, labels = plt.xticks(x_axis, x_ticks, rotation="vertical")

plt.bar(x_axis, y_follow_up, color="black", label="LSST + follow-up")
plt.ylabel("lensed SNe Ia with good delay")
plt.legend()
plt.tight_layout()
plt.savefig("%sMetric.png" % (plot_path), dpi=300)
