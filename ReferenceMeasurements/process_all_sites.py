#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    process_all_sites.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        03/12/2023 16:59

# this script automatically estimates DHRs using all reference measurements available at the input directory

from cal_diffuse_radiation import *
import pandas as pd
import os

if __name__ == "__main__":

    Filtered_RM1_DIR = './FilterRawData'
    RM1_DIR = './RM1_data'

    # Step1: For sites do not have diffuse radiation measurements available, make an estimation using model.
    # check if the file is csv. If yes, load variables and values from the csv
    for filename in os.listdir(RM1_DIR):
        if filename.endswith('.csv'):
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
            # get the latitude and longitude from the csv file data in same file name, but in directory of RM1_DIR
            raw_data = pd.read_csv(RM1_DIR + '/' + filename, delimiter=';')
            lat = raw_data['Lat_IS'][0]
            lon = raw_data['Lon_IS'][0]
            print('Estimating diffuse radiation for site: {} at lat: {} and lon: {}'.format(filename, lat, lon))


        else:
            continue



