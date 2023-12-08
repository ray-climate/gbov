#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    process_all_sites.py
@Author:      Dr. Rui Song
@Email:       rui.song@physics.ox.ac.uk
@Time:        03/12/2023 16:59

This script automatically estimates Diffuse Horizontal Radiation (DHR) using reference measurements available in the
input directory."""
import numpy as np

from cal_dhr_bhr import *
import pandas as pd
import datetime
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plot_dhr(dhr_values, dates, output_path):
    """
    Plots DHR values against measurement dates.

    :param dhr_values: List of DHR values.
    :param dates: List of measurement dates in string format '%Y-%m-%d'.
    :param output_path: Path to save the plot.
    """
    dates = [datetime.datetime.strptime(date, "%Y-%m-%d") for date in dates]

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.scatter(dates, dhr_values, label='DHR', color='red', marker='+', s=80)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45)

    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(16)

    ax.set_ylim([0, 0.5])
    ax.set_xlabel('Date', fontsize=16)
    ax.set_ylabel('DHR', fontsize=16)
    ax.legend(fontsize=16)
    ax.set_title(output_path.split('/')[-1][:-4], fontsize=16)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_bhr(dhr_values, dates, output_path):
    """
    Plots BHR values against measurement dates.

    :param bhr_values: List of DHR values.
    :param dates: List of measurement dates in string format '%Y-%m-%d'.
    :param output_path: Path to save the plot.
    """
    dates = [datetime.datetime.strptime(date, "%Y-%m-%d") for date in dates]

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.scatter(dates, dhr_values, label='BHR', color='red', marker='+', s=80)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45)

    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(16)

    ax.set_ylim([0, 0.5])
    ax.set_xlabel('Date', fontsize=16)
    ax.set_ylabel('BHR', fontsize=16)
    ax.legend(fontsize=16)
    ax.set_title(output_path.split('/')[-1][:-4], fontsize=16)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_tocr(dhr_values, dates, output_path):
    """
    Plots DHR values against measurement dates.

    :param dhr_values: List of DHR values.
    :param dates: List of measurement dates in string format '%Y-%m-%d'.
    :param output_path: Path to save the plot.
    """
    dates = [datetime.datetime.strptime(date, "%Y-%m-%d") for date in dates]

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.scatter(dates, dhr_values, label='ToC-R', color='red', marker='+', s=80)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45)

    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(16)

    ax.set_ylim([0, 0.5])
    ax.set_xlabel('Date', fontsize=16)
    ax.set_ylabel('ToC-R', fontsize=16)
    ax.legend(fontsize=16)
    ax.set_title(output_path.split('/')[-1][:-4], fontsize=16)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def main():
    """
    Main function to process all sites for DHR estimation.
    """
    tower_dir = '../sites'
    filtered_rm1_dir = './FilterRawData'
    rm1_dir = './RM1_data'
    save_dhr_dir = './OUTPUT_dhr_bhr_tocr'

    if not os.path.exists(save_dhr_dir):
        os.makedirs(save_dhr_dir)

    alpha_dhr = 0.15
    alpha_bhr = 0.85
    print('Start processing all sites'')
    for filename in os.listdir(rm1_dir):
        if filename.endswith('GBOV_RM01_Calperum_001_20211231T000000Z_20221231T140000Z_008_ACR_V2.0.csv'):

            dhr_values, bhr_values, albedo_dates = [], [], []
            tocr_values, tocr_dates = [], []

            print(f'Processing file: {filename}')
            data = pd.read_csv(f'{filtered_rm1_dir}/{filename}', delimiter=',')
            data['Datetime'] = pd.to_datetime(data['Datetime'])
            data[['SW_IN', 'SW_OUT', 'SW_DIF']] = data[['SW_IN', 'SW_OUT', 'SW_DIF']].apply(pd.to_numeric, errors='coerce')

            raw_data = pd.read_csv(f'{rm1_dir}/{filename}', delimiter=';')
            lat, lon = raw_data.at[0, 'Lat_IS'], raw_data.at[0, 'Lon_IS']
            data['sza'] = get_sza(data['Datetime'], lat, lon)
            data['SW_DIF_ratio'] = data['SW_DIF'] / data['SW_IN']

            tower_dates = np.loadtxt(f'{tower_dir}/central_dates.txt', dtype=object, delimiter=",", skiprows=1)
            start_year, end_year = filename.split('_')[4][:4], filename.split('_')[5][:4]
            year_range = range(int(start_year), int(end_year) + 1)
            csv_start_date, csv_end_date = filename.split('_')[4][:8], filename.split('_')[5][:8]

            for year in year_range:
                for m, (start_mdy, nominal_mdy, end_mdy) in enumerate(tower_dates):

                    start_date = f"{year - 1 if m <= 1 else year}-{start_mdy}"
                    end_date = f"{year + 1 if m == len(tower_dates) - 1 else year}-{end_mdy}"
                    nominal_date = f"{year}-{nominal_mdy}"

                    # check if nominal date is within the start and end date. nominal_date is YYYY-MM-DD, start and end date are YYYYMMDD
                    if not (int(csv_start_date.replace('-', '')) <= int(nominal_date.replace('-', '')) <= int(csv_end_date.replace('-', ''))):
                        continue

                    print(
                        f"Calculating albedo for {filename} using start date {start_date}, nominal date {nominal_date} "
                        f"and end date {end_date}")

                    try:
                        toc_r = cal_tower_tocr(data, start_date, nominal_date, end_date, lat, lon, alpha_dhr, alpha_bhr)

                        tocr_dates.append(nominal_date)
                        tocr_values.append(toc_r)

                    except:
                        print('NO enough data for TOC_R calculation')

                    try:
                        (dhr_value, bhr_value) = cal_tower_dhr_bhr(data, start_date, nominal_date, end_date, lat, lon, alpha_dhr, alpha_bhr)

                        albedo_dates.append(nominal_date)
                        dhr_values.append(dhr_value)
                        bhr_values.append(bhr_value)

                    except:
                        print('NO enough data for DHR BHR calculation')

            # if number of valid tocr_values is less than 10, then interpolate between the first and last valid tocr_values based on tocr_dates
            # tocr_dates is a list of strings, in the format of YYYY-MM-DD. Convert before using np.interp
            if len(np.asarray(tocr_values)[~np.isnan(tocr_values)]) < 12:

                print('Interpolating TOC_R values')
                tocr_dates_datetime = [datetime.datetime.strptime(date, "%Y-%m-%d") for date in tocr_dates]
                tocr_dates_datetime = mdates.date2num(tocr_dates_datetime)

                # Identify indices of first and last non-NaN values
                valid_indices = np.where(~np.isnan(tocr_values))[0]
                first_valid_index = valid_indices[0]
                last_valid_index = valid_indices[-1]
                print(f'First valid index: {first_valid_index}, last valid index: {last_valid_index}')
                print(tocr_dates[first_valid_index], tocr_dates[last_valid_index])
                # Perform interpolation only between the first and last valid values
                tocr_values[first_valid_index:last_valid_index + 1] = np.interp(
                    tocr_dates_datetime[first_valid_index:last_valid_index + 1],
                    tocr_dates_datetime[valid_indices],
                    np.asarray(tocr_values)[valid_indices])

            # if number of valid dhr_values is less than 10, then interpolate between the first and last valid dhr_values based on albedo_dates
            # albedo_dates is a list of strings, in the format of YYYY-MM-DD. Convert before using np.interp
            if len(np.asarray(dhr_values)[~np.isnan(dhr_values)]) < 12:

                    print('Interpolating DHR values')
                    albedo_dates_datetime = [datetime.datetime.strptime(date, "%Y-%m-%d") for date in albedo_dates]
                    albedo_dates_datetime = mdates.date2num(albedo_dates_datetime)

                    # Identify indices of first and last non-NaN values
                    valid_indices = np.where(~np.isnan(dhr_values))[0]
                    first_valid_index = valid_indices[0]
                    last_valid_index = valid_indices[-1]
                    print(f'First valid index: {first_valid_index}, last valid index: {last_valid_index}')
                    print(albedo_dates[first_valid_index], albedo_dates[last_valid_index])
                    # Perform interpolation only between the first and last valid values
                    dhr_values[first_valid_index:last_valid_index + 1] = np.interp(
                        albedo_dates_datetime[first_valid_index:last_valid_index + 1],
                        albedo_dates_datetime[valid_indices],
                        np.asarray(dhr_values)[valid_indices])

            # if number of valid bhr_values is less than 10, then interpolate between the first and last valid bhr_values based on albedo_dates
            # albedo_dates is a list of strings, in the format of YYYY-MM-DD. Convert before using np.interp
            if len(np.asarray(bhr_values)[~np.isnan(bhr_values)]) < 12:

                        print('Interpolating BHR values')
                        albedo_dates_datetime = [datetime.datetime.strptime(date, "%Y-%m-%d") for date in albedo_dates]
                        albedo_dates_datetime = mdates.date2num(albedo_dates_datetime)

                        # Identify indices of first and last non-NaN values
                        valid_indices = np.where(~np.isnan(bhr_values))[0]
                        first_valid_index = valid_indices[0]
                        last_valid_index = valid_indices[-1]
                        print(f'First valid index: {first_valid_index}, last valid index: {last_valid_index}')
                        print(albedo_dates[first_valid_index], albedo_dates[last_valid_index])
                        # Perform interpolation only between the first and last valid values
                        bhr_values[first_valid_index:last_valid_index + 1] = np.interp(
                            albedo_dates_datetime[first_valid_index:last_valid_index + 1],
                            albedo_dates_datetime[valid_indices],
                            np.asarray(bhr_values)[valid_indices])

            with open(f'{save_dhr_dir}/{filename[:-4]}_TOC_R.csv', 'w', newline='') as output:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerow(('Datetime', 'TOC_R'))
                for i in range(len(tocr_dates)):
                    writer.writerow((tocr_dates[i], tocr_values[i]))

            with open(f'{save_dhr_dir}/{filename[:-4]}_DHR.csv', 'w', newline='') as output:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerow(('Datetime', 'DHR'))
                for i in range(len(albedo_dates)):
                    writer.writerow((albedo_dates[i], dhr_values[i]))

            with open(f'{save_dhr_dir}/{filename[:-4]}_BHR.csv', 'w', newline='') as output:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerow(('Datetime', 'BHR'))
                for i in range(len(albedo_dates)):
                    writer.writerow((albedo_dates[i], bhr_values[i]))

            plot_dhr(dhr_values, albedo_dates, f'{save_dhr_dir}/{filename[:-4]}_DHR.png')
            plot_bhr(bhr_values, albedo_dates, f'{save_dhr_dir}/{filename[:-4]}_BHR.png')
            plot_tocr(tocr_values, tocr_dates, f'{save_dhr_dir}/{filename[:-4]}_TOC_R.png')


if __name__ == "__main__":
    main()
