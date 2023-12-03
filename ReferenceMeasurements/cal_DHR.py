#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_DHR.py
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

def cal_tower_dhr(data, start_datetime, nominal_datetime, end_datetime, latitude, longitude, alpha_DHR):

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
    data = data[data['SW_DIF_ratio'] < alpha_DHR]
    dhr_avg = data['SW_OUT'] / data['SW_IN']

    if len(dhr_avg) > 0:
        DHR_value = np.mean(dhr_avg[dhr_avg>0])
    else:
        DHR_value = np.nan
    print("DHR is calculated: %s" % DHR_value)

    return DHR_value