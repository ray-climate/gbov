#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    process_all_sites.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        03/12/2023 16:59

# this script automatically estimates DHRs using all reference measurements available at the input directory

from cal_diffuse_radiation import *
from cal_DHR import *
import pandas as pd
import datetime
import os

if __name__ == "__main__":

    TOWER_DIR = '../sites'
    Filtered_RM1_DIR = './FilterRawData'
    RM1_DIR = './RM1_data'
    alpha_DHR = 0.15

    albedo_date_all = []
    dhr_value_all = []

    # Step1: For sites do not have diffuse radiation measurements available, make an estimation using model.
    # check if the file is csv. If yes, load variables and values from the csv
    for filename in os.listdir(RM1_DIR):
        if filename.endswith('GBOV_RM01_Brasschaat_001_20210101T000000Z_20211231T233000Z_007_ACR_V2.0.csv'):
            print('Processing file: {}'.format(filename))
            # read the csv file
            data = pd.read_csv(Filtered_RM1_DIR + '/' + filename, delimiter=',')
            # Convert Datetime to datetime
            data['Datetime'] = pd.to_datetime(data['Datetime'], format='%Y-%m-%dT%H:%M:%S.%f')
            # Convert columns to numeric, replacing non-numeric values with NaN
            for col in ['SW_IN', 'SW_OUT', 'SW_DIF']:
                data[col] = pd.to_numeric(data[col], errors='coerce')
            # get the time, SW_IN and SW_OUT
            time = data['Datetime']
            SW_IN = data['SW_IN']
            SW_OUT = data['SW_OUT']
            SW_DIF = data['SW_DIF']

            # get the latitude and longitude from the csv file data in same file name, but in directory of RM1_DIR
            raw_data = pd.read_csv(RM1_DIR + '/' + filename, delimiter=';')
            lat = raw_data['Lat_IS'][0]
            lon = raw_data['Lon_IS'][0]

            # add a column of solar zenith angle to the data
            sza = get_sza(time, lat, lon)
            data['sza'] = sza

            # add a column of DIF_ratio to the data
            SW_DIF_ratio = SW_DIF / SW_IN
            data['SW_DIF_ratio'] = SW_DIF_ratio

            # get the central processing dates for processing
            tower_central_dates_data = np.loadtxt(TOWER_DIR + '/central_dates.txt', dtype=object, delimiter=",", skiprows=1)
            # get the starting and ending processing year from the filename of the csv file
            start_year = filename.split('_')[4][0:4]
            end_year = filename.split('_')[5][0:4]
            year_range = np.arange(int(start_year), int(end_year) + 1, 1)

            for year_j in year_range:
                for m in range(tower_central_dates_data.shape[0]):
                    start_date = str(year_j) + '-' + tower_central_dates_data[m, 0]
                    nominal_date = str(year_j) + '-' + tower_central_dates_data[m, 1]
                    end_date = str(year_j) + '-' + tower_central_dates_data[m, 2]

                    if m <= 1:
                        start_date = str(year_j - 1) + '-' + tower_central_dates_data[m, 0]

                    if m == tower_central_dates_data.shape[0] - 1:
                        end_date = str(year_j + 1) + '-' + tower_central_dates_data[m, 2]

                    start_date_datetime = datetime.datetime.strptime('%s' % start_date, "%Y-%m-%d")
                    nominal_date_datetime = datetime.datetime.strptime('%s' % nominal_date, "%Y-%m-%d")
                    end_date_datetime = datetime.datetime.strptime('%s' % end_date, "%Y-%m-%d")

                    print("------------------------------>")
                    print("Calculating albedo for %s using start date %s, nominal date %s and end date %s" % (filename, start_date, nominal_date, end_date))

                    DHR_value = cal_tower_dhr(data, start_date_datetime, nominal_date_datetime, end_date_datetime, lat, lon, alpha_DHR)

                    albedo_date_all.append(nominal_date)
                    dhr_value_all.append(DHR_value)

        else:
            continue



