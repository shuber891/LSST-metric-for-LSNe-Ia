#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:30:24 2019

@author: shuber
"""

import numpy as np
from astropy import units as u
import csv
import matplotlib.pyplot as plt


# To run this code you have to download the files from https://me.lsst.eu/gris/
# and the opsim_runs_strategies_input.csv file from 
# https://github.com/shuber891/LSST-metric-for-LSNe-Ia
# Further chance the data_path and output_path files.

# Choose the right version that it works 13 for v1.3 and 14 for v1.4
version = 14

# This is where your cadence strategy file is located
data_path = "/afs/mpa/data/shuber/data/LSST_MAF_data/SL_fbs%i/" % version
output_path = "/afs/mpa/home/shuber/research/MicrolensSN/plot/"


def f_get_median(array):
    l_median = np.ndarray.tolist(array)
    l_median.sort()
    median = l_median[int(len(l_median)*0.5)] * array.unit
    return median


def f_fit_function(gap_median_filter):
    # has been constructed by fitting 13 observing strategies
    fit = 2.15 * np.exp(0.37 * gap_median_filter.value)
    return fit

def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = np.remainder(RA+360-org,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection)
    ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
    plt.savefig("/afs/mpa/home/shuber/research/MicrolensSN/plot/%s.png" % title, dpi = 300)

observing_strategies = []
# to define which strategies to investigate
file_strategies = open("%sopsim_runs_strategies_input_%i.csv" % (data_path,version), "r")
filte_to_write_on = open("%sopsim_runs%i.csv" % (output_path,version), "w")
writer = csv.writer(filte_to_write_on)

counter_line = 0
for line in file_strategies:
    line = line.strip()
    if counter_line == 0:
        line_to_write = line.split(",")
        line_to_write.append("lensed SNe Ia with good time delay")
        writer.writerow(line_to_write)
    if counter_line >= 1:
        if version == 13:
            observing_strategy = line.split(",")[1][:-3]
        elif version == 14:
            observing_strategy = line.split(",")[1][:-7]
        #print observing_strategy

        array = np.load("%s%s_SL.npy" % (data_path, observing_strategy))

        gap_median_filter = f_get_median(array["gap_median"] * u.day)
        #cumulative_season_length = f_get_median((array["season_length"] * u.day).to(u.year)) 
        cumulative_season_length = np.mean((array["season_length"] * u.day).to(u.year))
        
        if observing_strategy == "rolling_mod2_sdf_0.10_v1.4_10yrs" or observing_strategy == "rolling_mod6_sdf_0.20_v1.4_10yrs" or observing_strategy == "footprint_stuck_rollingv1.4_10yrs" or observing_strategy == "baseline_v1.4_10yrs":
            """
            plt.figure(observing_strategy+"1")
            plt.scatter(dec,array["season_length"] / 365,label=observing_strategy)
            plt.xlabel("Dec [deg]")
            plt.ylabel("cumulative season length [yr]")
            plt.savefig("/afs/mpa/home/shuber/research/MicrolensSN/plot/csl_%s.png" % observing_strategy)
            plt.legend()
            
            
            plot_mwd(RA = ra, Dec = dec, title = observing_strategy) #ra in (0,360) deg, dec (-90,90)           
            """

        area_pixel = array["area"] * u.deg**2
        survey_area = np.sum(area_pixel)

        N_OM10 = 45.7
        survey_area_OM10 = 20000 * u.deg**2
        cumulative_season_length_OM10 = 2.5 * u.year

        # rescale numbers of OM10 with cumulative_season_length
        # and survey_area of a given cadence strategy
        total_number_of_LSNeIa_0M10 = N_OM10 * survey_area/survey_area_OM10 * cumulative_season_length/cumulative_season_length_OM10

        # this estimates the number of well LSNeIa which have a well
        # measured time delay, i.e. that at least one image as
        # accuracy < 1 percent and precision < 5 percent
        number_of_good_delay_LSNeIa = total_number_of_LSNeIa_0M10 / f_fit_function(gap_median_filter)

        #print observing_strategy, survey_area, cumulative_season_length, gap_median_filter, number_of_good_delay_LSNeIa

        line_to_write = line.split(",")
        line_to_write.append(number_of_good_delay_LSNeIa)
        writer.writerow(line_to_write)

    counter_line += 1

file_strategies.close()
filte_to_write_on.close()