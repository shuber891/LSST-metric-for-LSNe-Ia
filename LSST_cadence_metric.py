#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:30:24 2019

@author: shuber
"""

import numpy as np
from astropy import units as u
import csv

# To run this code you have to download the files from https://me.lsst.eu/gris/
# and the opsim_runs_strategies_input.csv file from ...
# Further chance the data_path and output_path files

# This is where your cadence strategy file is located
data_path = "/afs/mpa/data/shuber/data/LSST_MAF_data/"
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


observing_strategies = []
file_strategies = open("%sopsim_runs_strategies_input.csv" % data_path, "r")
filte_to_write_on = open("%sopsim_runs.csv" % output_path, "w")
writer = csv.writer(filte_to_write_on)

counter_line = 0
for line in file_strategies:
    line = line.strip()
    if counter_line == 0:
        line_to_write = line.split(",")
        line_to_write.append("Number of LSNe Ia")
        writer.writerow(line_to_write)
    if counter_line >= 1:
        print line.split(",")
        observing_strategy = line.split(",")[1][:-3]

        if observing_strategy == "baseline_v1.3_1yrs" or observing_strategy == "tde_illum75_v1.3_10yrs" or observing_strategy == "uer_illum75_v1.3_10yrs":
            line_to_write = line.split(",")
            line_to_write.append("-")
            writer.writerow(line_to_write)
        else:
            array = np.load("%s%s_SL.npy" % (data_path, observing_strategy))

            gap_median_filter = f_get_median(array["gap_median"] * u.day)
            cumulative_season_length = f_get_median((array["season_length"] * u.day).to(u.year)) 

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

            print observing_strategy, survey_area, cumulative_season_length, gap_median_filter, number_of_good_delay_LSNeIa

            line_to_write = line.split(",")
            line_to_write.append(number_of_good_delay_LSNeIa)
            writer.writerow(line_to_write)

    counter_line += 1

file_strategies.close()
filte_to_write_on.close()