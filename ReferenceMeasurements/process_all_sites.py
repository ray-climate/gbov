#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    process_all_sites.py
@Author:      Dr. Rui Song
@Email:       rui.song@physics.ox.ac.uk
@Time:        03/12/2023 16:59

This script automatically estimates Diffuse Horizontal Radiation (DHR) using reference measurements available in the
input directory."""

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

    ax.set_ylim([0, 0.3])
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

    ax.set_ylim([0, 0.3])
    ax.set_xlabel('Date', fontsize=16)
    ax.set_ylabel('BHR', fontsize=16)
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
    save_dhr_dir = './OUTPUT_dhr_bhr'

    if not os.path.exists(save_dhr_dir):
        os.makedirs(save_dhr_dir)

    alpha_dhr = 0.15
    alpha_bhr = 0.85

    for filename in os.listdir(rm1_dir):
        if filename.endswith('.csv'):

            dhr_values, bhr_values, albedo_dates = [], [], []

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

            for year in year_range:
                for m, (start_mdy, nominal_mdy, end_mdy) in enumerate(tower_dates):

                    start_date = f"{year - 1 if m <= 1 else year}-{start_mdy}"
                    end_date = f"{year + 1 if m == len(tower_dates) - 1 else year}-{end_mdy}"
                    nominal_date = f"{year}-{nominal_mdy}"

                    print(
                        f"Calculating albedo for {filename} using start date {start_date}, nominal date {nominal_date} "
                        f"and end date {end_date}")

                    try:
                        (dhr_value, bhr_value) = cal_tower_dhr_bhr(data, start_date, nominal_date, end_date, lat, lon, alpha_dhr, alpha_bhr)

                        albedo_dates.append(nominal_date)
                        dhr_values.append(dhr_value)
                        bhr_values.append(bhr_value)

                    except:
                        print('NO enough data for DHR BHR calculation')

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


if __name__ == "__main__":
    main()
