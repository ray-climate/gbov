#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_dhr_bhr.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        03/12/2023 21:24

from astral.location import Location
import numpy as np
import datetime
import os.path
import pvlib
import pytz
import csv
import os

def get_sza(dates, latitude, longitude):

    szaMax = 75.  # the applied BRDF model is only valid for the maximum sza of 75.
    szaGap = 5.

    # solar zenith angle. column 6.
    sza = np.zeros((dates.shape[0]))
    for k in range(len(dates)):
        sun_angles = pvlib.solarposition.get_solarposition(dates[k], latitude, longitude)
        sza[k] = sun_angles.zenith.values[0]

    sza = np.asarray(sza)
    sza = sza.reshape((sza.size, 1))

    return sza

def cal_tower_dhr_bhr(data, start_datetime, nominal_datetime, end_datetime, latitude, longitude, alpha_DHR, alpha_BHR):

    szaMax = 75.  # the applied BRDF model is only valid for the maximum sza of 75.
    szaGap = 5.

    # data filtering (date and time)
    central_time = np.copy(nominal_datetime)
    start_time = np.copy(start_datetime)
    end_time = np.copy(end_datetime)

    # drop data (in panda dataframe) with variable 'Datetime' larger than the start_time
    print('size of raw data: ', data.shape[0])
    data = data[(data['Datetime'] > start_datetime) & (data['Datetime'] <= end_datetime)]
    print('filter datetime', data.shape[0])

    if data.shape[0] == 0:
        return -9999

    # data filtering (solar zenith angle)
    data = data[data['sza'] < szaMax]
    print('filter sza', data.shape[0])

    # data filtering (valid SW_IN)
    data = data[data['SW_IN'] > 30]
    print('filter valid SW_IN', data.shape[0])

    # data filtering (valid SW_OUT)
    data = data[data['SW_OUT'] > 30]
    print('filter valid SW_OUT', data.shape[0])

    ###########################
    # extract DHR information #
    ###########################

    print("DHR is calculated by averaging data with diffuse ratio smaller than %s" % alpha_DHR)
    # cal dhr_avg using data['SW_IN'] / data['SW_OUT'], but only when data['SW_DIF_ratio'] < alpha_DHR

    while alpha_DHR <= 0.25:

        data_dhr = data[data['SW_DIF_ratio'] <= alpha_DHR]
        dhr_avg = data_dhr['SW_OUT'] / data_dhr['SW_IN']

        if len(dhr_avg[dhr_avg > 0]) > 0:
            DHR_value = np.mean(dhr_avg[dhr_avg > 0])
            break
        else:
            alpha_DHR += 0.02

    # set DHR_value to nan if alpha_DHR > 0.25
    if alpha_DHR > 0.25:
        DHR_value = np.nan

    print("DHR is calculated: %s" % DHR_value)

    while alpha_BHR >= 0.5:

        data_bhr = data[data['SW_DIF_ratio'] >= alpha_BHR]
        bhr_avg = data_bhr['SW_OUT'] / data_bhr['SW_IN']

        if len(bhr_avg[bhr_avg > 0]) > 0:
            BHR_value = np.mean(bhr_avg[bhr_avg>0])
            break
        else:
            alpha_BHR -= 0.02

    # set BHR_value to nan if alpha_BHR < 0.75
    if alpha_BHR < 0.5:
        BHR_value = np.nan

    print("BHR is calculated: %s" % BHR_value)

    return DHR_value, BHR_value