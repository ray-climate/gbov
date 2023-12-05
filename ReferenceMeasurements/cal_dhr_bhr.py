#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_dhr_bhr.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        03/12/2023 21:24

from timezonefinder import TimezoneFinder
from astral.location import Location
import pandas as pd
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

    # data filtering (valid SW_DIF_ratio)
    data = data[data['SW_DIF_ratio'] > 0]
    print('filter valid SW_DIF_ratio', data.shape[0])

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

def cal_tower_tocr(data, start_datetime, nominal_datetime, end_datetime, latitude, longitude, alpha_DHR, alpha_BHR):

    # look up table.
    I = {}
    I[0] = [1.0, -0.997910, -0.00894619]
    I[5] = [1.0, -0.998980, -0.00837790]
    I[10] = [1.0, -1.00197, -0.00665391]
    I[15] = [1.0, -1.00702, -0.00371872]
    I[20] = [1.0, -1.01438, 0.000524714]
    I[25] = [1.0, -1.02443, 0.00621877]
    I[30] = [1.0, -1.03773, 0.0135606]
    I[35] = [1.0, -1.05501, 0.0228129]
    I[40] = [1.0, -1.07742, 0.0343240]
    I[45] = [1.0, -1.10665, 0.0485505]
    I[50] = [1.0, -1.14526, 0.0661051]
    I[55] = [1.0, -1.19740, 0.0878086]
    I[60] = [1.0, -1.27008, 0.114795]
    I[65] = [1.0, -1.37595, 0.148698]
    I[70] = [1.0, -1.54059, 0.191944]
    I[75] = [1.0, -1.82419, 0.248471]
    I[80] = [1.0, -2.40820, 0.325351]
    I[85] = [1.0, -4.20369, 0.438371]
    I['bhr'] = [1.0, -1.28159, 0.0802838]

    szaMax = 75.  # the applied BRDF model is only valid for the maximum sza of 75.
    szaGap = 5.

    # BRDF function 1
    def brdf_f1(sza, vza, phi):

        sza = np.deg2rad(sza)
        vza = np.deg2rad(vza)
        phi = np.deg2rad(phi)

        parameter_1 = 1. / (2. * np.pi) * ((np.pi - phi) * np.cos(phi) + np.sin(phi)) * np.tan(sza) * np.tan(vza)
        parameter_2 = 1. / np.pi * (np.tan(sza) + np.tan(vza) + np.sqrt(
            np.tan(sza) ** 2 + np.tan(vza) ** 2 - 2 * np.tan(sza) * np.tan(vza) * np.cos(phi)))

        return (parameter_1 - parameter_2)

    # BRDF function 2
    def brdf_f2(sza, vza, phi):
        sza = np.deg2rad(sza)
        vza = np.deg2rad(vza)
        phi = np.deg2rad(phi)

        ci = np.arccos(np.cos(sza) * np.cos(vza) + np.sin(vza) * np.sin(sza) * np.cos(phi))

        return 4. / (3. * np.pi) / (np.cos(sza) + np.cos(vza)) * (
                    (np.pi / 2. - ci) * np.cos(ci) + np.sin(ci)) - 1. / 3.

    # convert local time to UTC, in column 1.
    tf = TimezoneFinder()
    tz = tf.timezone_at(lng=longitude, lat=latitude)
    tzInfo = pytz.timezone(tz)

    dates_local_format = pd.to_datetime(data['Datetime'])
    dates_local_format = [i.replace(tzinfo=tzInfo).astimezone(tzInfo) for i in dates_local_format]
    data['Datetime_local'] = dates_local_format

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

    # data filtering (valid SW_DIF_ratio)
    data = data[data['SW_DIF_ratio'] > 0]
    data = data[data['SW_DIF_ratio'] < 1.]
    print('filter valid SW_DIF_ratio', data.shape[0])

    ###########################
    # extract BHR information #
    ###########################

    BHR_value = np.nan
    BHR_beta = np.nan

    while alpha_BHR >= 0.7:

        data_bhr = data[data['SW_DIF_ratio'] >= alpha_BHR]
        bhr_avg = data_bhr['SW_OUT'] / data_bhr['SW_IN']

        if len(bhr_avg[bhr_avg > 0]) > 0:
            BHR_value = np.mean(bhr_avg[bhr_avg > 0])
            BHR_beta = np.mean(data_bhr['SW_DIF_ratio'])
            break
        else:
            alpha_BHR -= 0.02

    BlueSkyAlbedo = data['SW_OUT'] / data['SW_IN']

    # set BHR_value to nan if alpha_BHR < 0.75
    if (BHR_value > 0.) & (BHR_beta > 0.):
        data['DHR_corrected'] = (BlueSkyAlbedo - BHR_beta * BHR_value) / (1 - BHR_beta)
    else:
        print("Not enough measurements to correct DHR")
        data['DHR_corrected'] = BlueSkyAlbedo

    print('DHR_corrected', data['DHR_corrected'])

    #######################################################################################
    # compute BRDF kernels using DHRs within the time window.
    #######################################################################################

    szaArr = np.array(range(0, int(szaMax) + 5, 5))

    wangle = 2.5
    num = len(szaArr)
    Y = np.zeros(num)
    X = np.zeros((szaArr.shape[0], 3))

    albedo = np.asarray(data['DHR_corrected'])
    zen = np.asarray(data['sza'])

    szaArr2 = []
    albedo_avg = {}
    i = 0

    for a in szaArr:
        w = np.where(np.abs(zen - a) <= wangle)
        if len(w[0]) > 0:
            szaArr2.append(a)
            if len(w[0]) > 1:
                Y[i] = np.mean(albedo[w])
                albedo_avg[a] = np.mean(albedo[w])
            else:
                Y[i] = albedo[w]
                albedo_avg[a] = albedo[w]
            X[i, 0], X[i, 1], X[i, 2] = I[a][0], I[a][1], I[a][2]
        i += 1

    try:
        Z = np.dot(X.T, X)
        Z = np.linalg.inv(Z)
        K = np.dot(Z, X.T)
        K = np.dot(K, Y)
        print("X is:", X)
        print("Y is:", Y)
        print("K is:", K)
        print("------------------------------>")
    except:
        print("Not enough measurements to calculate BRDF kernels")
        return np.nan

    #######################################################################################
    # Compute TOC-R at specified time from BRDF Kernels
    #######################################################################################

    central_time_as_datetime = datetime.datetime.strptime(str(central_time), '%Y-%m-%d')
    tocr_time = central_time_as_datetime + datetime.timedelta(hours=10)
    tocr_time_local_format = tzInfo.localize(tocr_time)
    tocr_time_UTC_format = tocr_time_local_format.astimezone(pytz.utc)

    print('ToC-R is calculated at local time: %s' % tocr_time_local_format)
    print('ToC-R is calculated at UTC time: %s' % tocr_time_UTC_format)

    tocr_angles = pvlib.solarposition.get_solarposition(tocr_time_UTC_format, latitude, longitude)
    sza_tocr = tocr_angles.zenith.values[0]

    print("solar zenith angle at 10:00 is %s " % sza_tocr)
    F = np.array([1, brdf_f1(sza_tocr, 0, 0), brdf_f2(sza_tocr, 0, 0)], 'f4')
    tocr_value = np.dot(K, F)

    if (tocr_value < 0) | (tocr_value > .5):
        tocr_value = np.nan

    print("ToC-R is calculated: %s" % tocr_value)

    return tocr_value