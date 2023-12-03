#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_diffuse_radiation.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        03/12/2023 17:09

import numpy as np

def cal_diff_radiation(lat, lon, time, SW_IN):
    """
    This function estimates the diffuse radiation if they are
    missing from the tower measurements.

    :param lat: latitude of the site.
    :param lon: longitude of the site, +/- on East/West hemisphere.
    :param time: local time of the measurement, format: datetime string.
    :param SW_IN: incoming shortwave solar radiation.

    :return:

    """

    if len(time) == 1:
        local_hour = time.hour
        local_minute = time.minute
        local_second = time.second
        local_standard_time = time.time
        local_number_of_day = time.timetuple().tm_yday
    else:
        local_hour = [local_time_i.hour for local_time_i in time]
        local_minute = [local_time_i.minute for local_time_i in time]
        local_second = [local_time_i.second for local_time_i in time]
        local_standard_time = [local_time_i.time() for local_time_i in time]
        local_number_of_day = [local_time_i.timetuple().tm_yday for local_time_i in time]

    local_hour = np.asarray(local_hour)
    local_minute = np.asarray(local_minute)
    local_second = np.asarray(local_second)
    local_number_of_day = np.asarray(local_number_of_day)

    # latitude is represented by varphi in radians.
    varphi = np.copy(np.radians(lat))

    # calculate declination angle (delta) in radians. Maleki et al. (2017), Section 2.
    Gamma = 2 * np.pi * (local_number_of_day - 1.) / 365.  # Gamma is the day angle in radians.
    delta = np.radians(23.45 * np.sin(np.radians(360. * (284. + local_number_of_day) / 365.)))

    # calculate hour angle (omega).
    B = 360. * (local_number_of_day - 81.) / 365.
    # ET is equation of time given by Tasdemiroglu (1998).
    ET = 9.87 * np.sin(2 * B) - 7.53 * np.cos(B) - 1.5 * np.cos(B)
    # LSMT is local standard time meridian, in degrees.
    LSTM = 15 * round(lon / 15)
    # LT is local time converted into hours.
    LT = local_hour + local_minute / 60. + local_second / 3600.

    # TC is     time correction factor to correct for longitude position
    # within the actual time zone.
    TC = 4. * (lon - LSTM) + ET
    # ST is local solar time.
    ST = LT + TC / 60.
    # omega is hour angle.
    omega = np.radians(15. * (ST - 12))

    # calculate equation of time, Spencer (1971).
    # E0 gives the equation of time (in minutes).
    E0 = 1.000110 + 0.034221 * np.cos(Gamma) + 0.001280 \
         * np.sin(Gamma) + 0.000719 * np.cos(2 * Gamma) + \
         0.000077 * np.sin(2 * Gamma)

    # Gon is incident irradiation calculated for the day number.
    Gon = 1367. * E0

    # calculate solar zenith angle cos_theta
    cos_theta = np.cos(varphi) * np.cos(delta) * np.cos(omega) \
                + np.sin(varphi) * np.sin(delta)
    cos_theta[cos_theta < 0] = 0

    sza = np.arccos(cos_theta)

    # calculate extraterrestrial radiation on surface tangent.
    G_0 = Gon * cos_theta

    # calculate diffuse ratio.
    beta = (G_0 - SW_IN) / G_0

    # estimate diffuse radiation.
    SW_DIF = SW_IN * beta

    return SW_DIF