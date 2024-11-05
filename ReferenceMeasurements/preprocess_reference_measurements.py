#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    preprocess_reference_measurements.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        28/11/2023 18:02

import pandas as pd
import numpy as np
import datetime
import csv
import os
from IPython import embed

# write a function to plot SW_IN, SW_OUT and SW_DIF
def plot_SW(SW_IN, SW_OUT, SW_DIF, MeasureTime, filename):

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    fig, ax = plt.subplots(figsize=(20, 8))
    plt.scatter(MeasureTime, SW_IN, label='SW_IN', color='blue', s=5, alpha=0.3)
    plt.scatter(MeasureTime, SW_DIF, label='SW_OUT', color='red', s=5, alpha=0.3)
    plt.scatter(MeasureTime, SW_OUT, label='SW_DIF', color='green', s=5, alpha=0.3)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(16)
    ax.set_ylim([0, 1500])
    ax.set_xlabel('Date', fontsize=16)
    ax.set_ylabel('SW (W/m2)', fontsize=16)
    ax.legend(fontsize=16)
    # set title by the string in filename betweewn the second and third '_'
    title = filename.split('_')[2]
    ax.set_title(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def preprocess_raw(DATA_DIR, filename, SAVE_DIR, file):

    data = pd.read_csv(DATA_DIR + '/' + filename, delimiter=';')

    # Convert TIME_IS to datetime
    data['TIME_IS'] = pd.to_datetime(data['TIME_IS'], format='%Y%m%dT%H%M%SZ')

    # Convert columns to numeric, replacing non-numeric values with NaN
    for col in ['QC_SW_IN', 'QC_SW_OUT', 'QC_SW_DIF']:
        data[col] = pd.to_numeric(data[col], errors='coerce')
        data[col].replace(1.0, np.nan, inplace=True)
    for col in ['SW_IN', 'SW_OUT', 'SW_DIF']:
        data[col] = pd.to_numeric(data[col], errors='coerce')
        data[col].replace(-999., np.nan, inplace=True)

    # Drop the rows with missing data
    data.dropna(subset=["SW_IN", "SW_OUT", "SW_DIF", "QC_SW_IN", "QC_SW_OUT", "QC_SW_DIF"], inplace=True)

    MeasureTime = np.copy(data['TIME_IS'])
    SW_IN = np.asarray(data['SW_IN'])
    SW_OUT = np.asarray(data['SW_OUT'])
    SW_DIF = np.asarray(data['SW_DIF'])

    SW_IN[SW_IN < 30.] = np.nan
    SW_OUT[SW_OUT < 30.] = np.nan
    SW_DIF[SW_DIF < 30.] = np.nan

    plot_SW(SW_IN, SW_OUT, SW_DIF, MeasureTime, SAVE_DIR + '/%s.png' % file[:-4])

    with open(SAVE_DIR + '/%s.csv' % filename[:-4], 'w') as output:

        writer = csv.writer(output, lineterminator='\n')
        writer.writerow(('Datetime', 'SW_IN', 'SW_DIF', 'SW_OUT'))
        for i in range(len(MeasureTime)):
            writer.writerow((MeasureTime[i], SW_IN[i], SW_DIF[i], SW_OUT[i]))

def main():
    RM1_DATA_DIR = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/GBOV-RM1-LP/RM01'
    SAVE_DIR = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/GBOV-RM1-LP/FilterRawData'

    filtered_filelist = []
    for file in os.listdir(SAVE_DIR):
        if file.endswith('csv'):
            filtered_filelist.append("_".join(file.split("_")[:9]))

    # Create output directory if not exists
    os.makedirs(SAVE_DIR, exist_ok=True)
    for file in os.listdir(RM1_DATA_DIR):
        not_yet_filtered = True
        for filtered_file in filtered_filelist:
            if filtered_file.split('_')[2] == file.split('_')[2]: not_yet_filtered = False

        if file.endswith('csv') and not_yet_filtered:
            print('Preprocessing file: {}'.format(file))
            preprocess_raw(RM1_DATA_DIR, file, SAVE_DIR, file)
            print('Preprocessing finished, saved to: {}'.format(SAVE_DIR + '/%s.csv' % file[:-4]))

if __name__ == "__main__":
    main()