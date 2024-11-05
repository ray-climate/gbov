#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_upscaling.py
# @Author:      Dr. Jean-Sébastien Carrière
# @Email:       jean-sebastien.carriere@acri-st.fr
# @Time:        10/10/2024 13:48

import os
import csv
import glob
import numpy as np
import pandas as pd
from osgeo import gdal
from random import random
from IPython import embed
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from collections import defaultdict
plt.rcParams.update({'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsfonts}'})


class Colors:
    RESET = '\033[0m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def plot_output(lastUpscale_dir, Xupscale, Yupscale, Yupscale_unc, Xgbov, Ygbov, Ygbov_unc, LEGEND, SAVE_LABEL):
    if np.isnan(np.array(Yupscale).astype(float)).all():
        print(Colors.RED + SAVE_LABEL + " upscaling failed" + Colors.RESET)
    elif not os.path.isfile(os.path.join(lastUpscale_dir, f"QCheck/Macro_{SAVE_LABEL}.png")):
        Ygbov = np.where(np.array(Ygbov) == "", np.nan, np.array(Ygbov))
        Ygbov_unc = np.where(np.array(Ygbov_unc) == "", np.nan, np.array(Ygbov_unc))
        startTime = datetime(2018, 1, 1, 0, 0)
        stopTime = datetime(2023, 12, 31, 0, 0)

        Ygbov, Yupscale = np.array(Ygbov).astype(float), np.array(Yupscale).astype(float)
        Ygbov_unc, Yupscale_unc = np.array(Ygbov_unc).astype(float), np.array(Yupscale_unc).astype(float)
        Xupscale = np.array(Xupscale)
        SortUpscale = Xupscale.argsort()
        Yupscale = Yupscale[SortUpscale, :]
        Yupscale_unc = Yupscale_unc[SortUpscale, :]
        Xupscale = Xupscale[SortUpscale]

        #################################
        ####### Nine plots
        #################################

        max_Y = np.nanmax(np.array(Yupscale).astype(float) + np.array(Yupscale_unc).astype(float))
        min_Y = np.nanmin(np.array(Yupscale).astype(float) - np.array(Yupscale_unc).astype(float))

        fig, ax = plt.subplots(nrows=9, figsize=(25, 45), dpi=200)

        for i in range(9):
            ax[i].spines['bottom'].set_color('#AAAAAA')
            ax[i].spines['top'].set_color('#AAAAAA')
            ax[i].spines['right'].set_color('#AAAAAA')
            ax[i].spines['left'].set_color('#AAAAAA')
            ax[i].tick_params(axis='x', colors='#444444')
            ax[i].tick_params(axis='y', colors='#444444')
            ax[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            ax[i].xaxis.set_major_locator(mdates.MonthLocator(interval=3))  # Every 3 months

            ax[i].errorbar(Xupscale, Yupscale[:, i], yerr=Yupscale_unc[:, i], color="black", fmt='o', capsize=5)

            ax[i].set_xlim(datetime(2023, 1, 1, 0, 0),
                        datetime(2023, 12, 31, 0, 0))
            ax[i].set_ylim(min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y))
            ax[i].set_xlabel(r"${\rm Time}$")
            ax[i].set_ylabel(r"${\rm " + LEGEND + "}$")
            ax[i].set_title(r"${\rm " + LEGEND + "}$")
            ax[i].grid(True, alpha=0.5)
            ax[i].set_facecolor("white")

            ax[i].plot([datetime(2023, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                    [0, 0], linestyle="--", alpha=0.5, color="black")
            ax[i].plot([datetime(2023, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                    [1, 1], linestyle="--", alpha=0.5, color="black")

            # Make every vertical gridline for January 1st bold
            for label in ax[i].get_xticklabels():
                date_str = label.get_text().replace('{\\fontfamily{\\familydefault}\\selectfont ', '').replace('}', '')
                date = datetime.strptime(date_str, "%Y-%m-%d")
                if date.month == 1 and date.day == 1:
                    x_pos = mdates.date2num(date)
                    ax[i].axvline(x=x_pos, color='orange', linewidth=2)

        fig.savefig(os.path.join(lastUpscale_dir, f"QCheck/{SAVE_LABEL}_individuals"), bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

        #################################
        ####### Mean plots with mean unc
        #################################

        Ygbov = np.nanmean(Ygbov, axis=1)
        Yupscale = np.nanmean(Yupscale, axis=1)

        Ygbov_unc = np.sqrt(np.nansum(Ygbov_unc ** 2, axis=1))/9
        Yupscale_unc = np.sqrt(np.nansum(Yupscale_unc ** 2, axis=1))/9

        if np.isnan(Ygbov.astype(float)).all():
            max_Y = np.nanmax(np.array(Yupscale).astype(float)+np.array(Yupscale_unc).astype(float))
            min_Y = np.nanmin(np.array(Yupscale).astype(float)-np.array(Yupscale_unc).astype(float))
        else:
            max_Y = np.max([np.nanmax(np.array(Yupscale).astype(float)+np.array(Yupscale_unc).astype(float)), np.nanmax(Ygbov.astype(float)+np.array(Ygbov_unc).astype(float))])
            min_Y = np.min([np.nanmin(np.array(Yupscale).astype(float)-np.array(Yupscale_unc).astype(float)), np.nanmin(Ygbov.astype(float)-np.array(Ygbov_unc).astype(float))])

        # Step 1: Split Xgbov and Ygbov by year
        grouped_Xgbov = defaultdict(list)
        grouped_Ygbov = defaultdict(list)
        grouped_Ygbov_unc = defaultdict(list)

        # Iterate through Xgbov and group the corresponding Ygbov rows by year
        for i, date in enumerate(Xgbov):
            year = date.year
            grouped_Xgbov[year].append(date)
            grouped_Ygbov[year].append(Ygbov[i])
            grouped_Ygbov_unc[year].append(Ygbov_unc[i])

        # Convert lists in grouped_Ygbov from lists to numpy arrays
        for year in grouped_Ygbov:
            grouped_Ygbov[year] = np.array(grouped_Ygbov[year])
        # Convert lists in grouped_Ygbov from lists to numpy arrays
        for year in grouped_Ygbov_unc:
            grouped_Ygbov_unc[year] = np.array(grouped_Ygbov_unc[year])

        # Convert dictionary to lists of blocks for Xgbov and Ygbov
        Xgbov_blocks = list(grouped_Xgbov.values())
        Ygbov_blocks = list(grouped_Ygbov.values())
        Ygbov_unc_blocks = list(grouped_Ygbov_unc.values())

        fig, ax = plt.subplots(figsize=(25, 5), dpi=200)

        ax.spines['bottom'].set_color('#AAAAAA')
        ax.spines['top'].set_color('#AAAAAA')
        ax.spines['right'].set_color('#AAAAAA')
        ax.spines['left'].set_color('#AAAAAA')
        ax.tick_params(axis='x', colors='#444444')
        ax.tick_params(axis='y', colors='#444444')
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))  # Every 3 months
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        # Plot each block (each year's data) separately to avoid connecting across years
        for i in range(len(Xgbov_blocks)):
            ax.plot(Xgbov_blocks[i], Ygbov_blocks[i], color="black")
            ax.fill_between(Xgbov_blocks[i],
                            Ygbov_blocks[i] + Ygbov_unc_blocks[i],
                            Ygbov_blocks[i] - Ygbov_unc_blocks[i],
                            alpha=0.3, color="black")

        ax.plot(Xupscale, Yupscale, color="black")
        ax.fill_between(Xupscale,
                        Yupscale + Yupscale_unc,
                        Yupscale - Yupscale_unc,
                        alpha=0.3, color="black")

        ax.fill_betweenx([min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y)],
                         datetime(2023, 1, 1, 0, 0),
                         datetime(2023, 12, 31, 0, 0),
                         alpha=0.3, color="orange")

        ax.plot([datetime(2018, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                [0, 0], linestyle="--", alpha=0.5, color="black")
        ax.plot([datetime(2018, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                [1, 1], linestyle="--", alpha=0.5, color="black")

        ax.set_xlim(startTime, stopTime)
        ax.set_ylim(min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y))
        ax.set_xlabel(r"${\rm Time}$")
        ax.set_ylabel(r"${\rm " + LEGEND + "}$")
        ax.set_title(r"${\rm " + LEGEND + "}$")
        ax.grid(True, alpha=0.5)
        ax.set_facecolor("white")

        # Make every vertical gridline for January 1st bold
        for label in ax.get_xticklabels():
            date_str = label.get_text().replace('{\\fontfamily{\\familydefault}\\selectfont ', '').replace('}', '')
            date = datetime.strptime(date_str, "%Y-%m-%d")
            if date.month == 1 and date.day == 1:
                x_pos = mdates.date2num(date)
                ax.axvline(x=x_pos, color='orange', linewidth=2)

        fig.savefig(os.path.join(lastUpscale_dir, f"QCheck/Macro_{SAVE_LABEL}"), bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

def plot_QC(lastUpscale_dir, Xupscale, Yupscale, Xgbov, Ygbov, LEGEND, SAVE_LABEL):
    if np.isnan(np.array(Yupscale).astype(float)).all():
        print(Colors.RED + SAVE_LABEL + " upscaling failed" + Colors.RESET)
    elif not os.path.isfile(os.path.join(lastUpscale_dir, f"QCheck/{SAVE_LABEL}.png")):
        Ygbov = np.where(np.array(Ygbov) == "", np.nan, np.array(Ygbov))
        if np.isnan(Ygbov.astype(float)).all():
            max_Y = np.nanmax(np.array(Yupscale).astype(float))
            min_Y = np.nanmin(np.array(Yupscale).astype(float))
        else:
            max_Y = np.max([np.nanmax(np.array(Yupscale).astype(float)), np.nanmax(Ygbov.astype(float))])
            min_Y = np.min([np.nanmin(np.array(Yupscale).astype(float)), np.nanmin(Ygbov.astype(float))])
        startTime = datetime(2018, 1, 1, 0, 0)
        stopTime = datetime(2023, 12, 31, 0, 0)

        #################################
        ####### The nine points together
        #################################

        fig, ax = plt.subplots(figsize=(25, 5), dpi=200)

        ax.spines['bottom'].set_color('#AAAAAA')
        ax.spines['top'].set_color('#AAAAAA')
        ax.spines['right'].set_color('#AAAAAA')
        ax.spines['left'].set_color('#AAAAAA')
        ax.tick_params(axis='x', colors='#444444')
        ax.tick_params(axis='y', colors='#444444')
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))  # Every 3 months
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        for pixel in range(9):
            ax.scatter(Xgbov, Ygbov[:, pixel].astype(float), color="black", marker="o")
            ax.scatter(Xupscale, np.array(Yupscale)[:, pixel].astype(float), color="black", marker="o")

        ax.fill_betweenx([min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y)],
                         datetime(2023, 1, 1, 0, 0),
                         datetime(2023, 12, 31, 0, 0),
                         alpha=0.3, color="orange")

        ax.set_xlim(startTime, stopTime)
        ax.set_ylim(min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y))
        ax.set_xlabel(r"${\rm Time}$")
        ax.set_ylabel(r"${\rm "+LEGEND+"}$")
        ax.set_title(r"${\rm "+LEGEND+"}$")
        ax.grid(True, alpha=0.5)
        ax.set_facecolor("white")

        ax.plot([datetime(2018, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                [0, 0], linestyle="--", alpha=0.5, color="black")
        ax.plot([datetime(2018, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                [1, 1], linestyle="--", alpha=0.5, color="black")


        # Make every vertical gridline for January 1st bold
        for label in ax.get_xticklabels():
            date_str = label.get_text().replace('{\\fontfamily{\\familydefault}\\selectfont ', '').replace('}', '')
            date = datetime.strptime(date_str, "%Y-%m-%d")
            if date.month == 1 and date.day == 1:
                x_pos = mdates.date2num(date)
                ax.axvline(x=x_pos, color='orange', linewidth=2)

        fig.savefig(os.path.join(lastUpscale_dir, f"QCheck/{SAVE_LABEL}"), bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

        #################################
        ####### The mean with the 1-sigma
        #################################

        Ygbov, Yupscale = np.array(Ygbov).astype(float), np.array(Yupscale).astype(float)
        Xupscale = np.array(Xupscale)
        SortUpscale = Xupscale.argsort()
        Yupscale = Yupscale[SortUpscale, :]
        Xupscale = Xupscale[SortUpscale]

        # Step 1: Split Xgbov and Ygbov by year
        grouped_Xgbov = defaultdict(list)
        grouped_Ygbov = defaultdict(list)

        # Iterate through Xgbov and group the corresponding Ygbov rows by year
        for i, date in enumerate(Xgbov):
            year = date.year
            grouped_Xgbov[year].append(date)
            grouped_Ygbov[year].append(Ygbov[i])

        # Convert lists in grouped_Ygbov from lists to numpy arrays
        for year in grouped_Ygbov:
            grouped_Ygbov[year] = np.array(grouped_Ygbov[year])

        # Convert dictionary to lists of blocks for Xgbov and Ygbov
        Xgbov_blocks = list(grouped_Xgbov.values())
        Ygbov_blocks = list(grouped_Ygbov.values())

        fig, ax = plt.subplots(figsize=(25, 5), dpi=200)

        ax.spines['bottom'].set_color('#AAAAAA')
        ax.spines['top'].set_color('#AAAAAA')
        ax.spines['right'].set_color('#AAAAAA')
        ax.spines['left'].set_color('#AAAAAA')
        ax.tick_params(axis='x', colors='#444444')
        ax.tick_params(axis='y', colors='#444444')
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))  # Every 3 months
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        # Plot each block (each year's data) separately to avoid connecting across years
        for i in range(len(Xgbov_blocks)):
            ax.errorbar(Xgbov_blocks[i], np.nanmean(Ygbov_blocks[i], axis=1), yerr=np.nanstd(Ygbov_blocks[i], axis=1), color="black", fmt='o', capsize=5)
            """ax.plot(Xgbov_blocks[i], np.nanmean(Ygbov_blocks[i], axis=1), color="black")
            ax.fill_between(Xgbov_blocks[i],
                            np.nanmean(Ygbov_blocks[i], axis=1) - np.nanstd(Ygbov_blocks[i], axis=1),
                            np.nanmean(Ygbov_blocks[i], axis=1) + np.nanstd(Ygbov_blocks[i], axis=1),
                            alpha=0.3, color="black")"""

        ax.errorbar(Xupscale, np.nanmean(Yupscale, axis=1), yerr=np.nanstd(Yupscale, axis=1), color="black", fmt='o', capsize=5)
        """ax.plot(Xupscale, np.nanmean(Yupscale, axis=1), color="black")
        ax.fill_between(Xupscale,
                        np.nanmean(Yupscale, axis=1) - np.nanstd(Yupscale, axis=1),
                        np.nanmean(Yupscale, axis=1) + np.nanstd(Yupscale, axis=1),
                        alpha=0.3, color="black")"""

        ax.fill_betweenx([min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y)],
                         datetime(2023, 1, 1, 0, 0),
                         datetime(2023, 12, 31, 0, 0),
                         alpha=0.3, color="orange")

        ax.plot([datetime(2018, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                [0, 0], linestyle="--", alpha=0.5, color="black")
        ax.plot([datetime(2018, 1, 1, 0, 0), datetime(2023, 12, 31, 0, 0)],
                [1, 1], linestyle="--", alpha=0.5, color="black")

        ax.set_xlim(startTime, stopTime)
        ax.set_ylim(min_Y - 0.1 * (max_Y - min_Y), max_Y + 0.1 * (max_Y - min_Y))
        ax.set_xlabel(r"${\rm Time}$")
        ax.set_ylabel(r"${\rm " + LEGEND + "}$")
        ax.set_title(r"${\rm " + LEGEND + "}$")
        ax.grid(True, alpha=0.5)
        ax.set_facecolor("white")

        # Make every vertical gridline for January 1st bold
        for label in ax.get_xticklabels():
            date_str = label.get_text().replace('{\\fontfamily{\\familydefault}\\selectfont ', '').replace('}', '')
            date = datetime.strptime(date_str, "%Y-%m-%d")
            if date.month == 1 and date.day == 1:
                x_pos = mdates.date2num(date)
                ax.axvline(x=x_pos, color='orange', linewidth=2)

        fig.savefig(os.path.join(lastUpscale_dir, f"QCheck/Macro_{SAVE_LABEL}"), bbox_inches='tight', pad_inches=0.1)
        plt.close('all')


        print(Colors.GREEN + SAVE_LABEL + " processed & saved" + Colors.RESET)
    else:
        print(Colors.GREEN + SAVE_LABEL + " already processed" + Colors.RESET)

lastUpscale_dir = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/GBOV-RM1-LP/LP01-02'
gbov_dir = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/GBOV-RM1-LP/GBOV_truth/'

embed()

Sites_list = [os.path.basename(f) for f in glob.glob(lastUpscale_dir+"/*") if "QCheck" not in f]
lon, lat = [], []
for SITE in Sites_list:
    print(Colors.BLUE  + SITE + Colors.RESET)
    #################################
    ####### Recover last upscale data
    #################################
    lastUpscale_files = glob.glob(lastUpscale_dir+f"/{SITE}/*.csv")
    Date_albedo_upscale, Date_TOCR_upscale = [], []
    TOCRupscale, TOCRUNCupscale = [], []
    BHRupscale, BHRUNCupscale, DHRupscale, DHRUNCupscale = [], [], [], []
    for lastUpscale in lastUpscale_files:
        if lastUpscale.split("_")[1] == "LP01":
            Date_TOCR_upscale.append(lastUpscale.split("_")[4])
            # Open and read a CSV file
            df = pd.read_csv(lastUpscale)
            TOCRupscale.append(df['TOC_R'].values)
            TOCRUNCupscale.append(df['TOC_R_unc'].values)

        elif lastUpscale.split("_")[1] == "LP02":
            Date_albedo_upscale.append(lastUpscale.split("_")[4])
            # Open and read a CSV file
            df = pd.read_csv(lastUpscale)
            BHRupscale.append(df['BHR'].values)
            DHRupscale.append(df['DHR'].values)
            BHRUNCupscale.append(df['BHR_unc'].values)
            DHRUNCupscale.append(df['DHR_unc'].values)
    Date_TOCR_upscale = [datetime.strptime(date, '%Y%m%d') for date in Date_TOCR_upscale]
    Date_albedo_upscale = [datetime.strptime(date, '%Y%m%d') for date in Date_albedo_upscale]

    lon.append(df['Longitude'].values)
    lat.append(df['Latitude'].values)

    #########################
    ####### Recover GBOV data
    #########################
    Date_TOCR_gbov, TOCRgbov, TOCRUNCgbov = [], [], []
    df = pd.read_csv(glob.glob(os.path.join(gbov_dir, "LP01", f'*macro*{SITE}*.csv'))[0])
    if SITE == "SouthernGreatPlains":
        initial_colum = df.columns[0]
        df = df.reset_index()
        df[initial_colum] = df['level_0']+df['level_1']+df[initial_colum]
        df = df[[initial_colum]]
    if SITE == "Hainich":
        initial_colum = df.columns[0]
        df = df.reset_index()
        df[initial_colum] = df['index']+df[initial_colum]
        df = df[[initial_colum]]
    df_tocr = df[df.columns[0]].str.split(';', expand=True)
    df_tocr.columns = df.columns[0].split(';')
    unique_dates = df_tocr['TIME_IS'].unique()
    for date_gbov in unique_dates:
        Date_TOCR_gbov.append(date_gbov.split('T')[0])
        TOCRgbov.append(df_tocr[df_tocr['TIME_IS'] == date_gbov]['ToC-R'].values)
        TOCRUNCgbov.append(df_tocr[df_tocr['TIME_IS'] == date_gbov]['ToC-R_unc'].values)
    Date_TOCR_gbov = [datetime.strptime(date, '%Y%m%d') for date in Date_TOCR_gbov]

    Date_ALBEDO_gbov, BHRgbov, BHRUNCgbov, DHRgbov, DHRUNCgbov = [], [], [], [], []
    df = pd.read_csv(glob.glob(os.path.join(gbov_dir, "LP02", f'*macro*{SITE}*.csv'))[0])
    if SITE == "SouthernGreatPlains":
        initial_colum = df.columns[0]
        df = df.reset_index()
        df[initial_colum] = df['level_0']+df['level_1']+df[initial_colum]
        df = df[[initial_colum]]
    if SITE == "Hainich":
        initial_colum = df.columns[0]
        df = df.reset_index()
        df[initial_colum] = df['index']+df[initial_colum]
        df = df[[initial_colum]]
    df_albedo = df[df.columns[0]].str.split(';', expand=True)
    df_albedo.columns = df.columns[0].split(';')
    unique_dates = df_albedo['TIME_IS'].unique()
    for date_gbov in unique_dates:
        Date_ALBEDO_gbov.append(date_gbov.split('T')[0])
        BHRgbov.append(df_albedo[df_albedo['TIME_IS'] == date_gbov]['BHR'].values)
        BHRUNCgbov.append(df_albedo[df_albedo['TIME_IS'] == date_gbov]['BHR_unc'].values)
        DHRgbov.append(df_albedo[df_albedo['TIME_IS'] == date_gbov]['DHR'].values)
        DHRUNCgbov.append(df_albedo[df_albedo['TIME_IS'] == date_gbov]['DHR_unc'].values)
    Date_ALBEDO_gbov = [datetime.strptime(date, '%Y%m%d') for date in Date_ALBEDO_gbov]

    #######################
    ####### Plot everything
    #######################
    plot_QC(lastUpscale_dir, Date_TOCR_upscale, TOCRupscale, Date_TOCR_gbov, TOCRgbov, "Top~of~canopy~reflectance", f"{SITE}_TOC-R")
    plot_QC(lastUpscale_dir, Date_TOCR_upscale, TOCRUNCupscale, Date_TOCR_gbov, TOCRUNCgbov, "Top~of~canopy~reflectance~uncertainty", f"{SITE}_TOC-R_unc")
    plot_QC(lastUpscale_dir, Date_albedo_upscale, BHRupscale, Date_ALBEDO_gbov, BHRgbov, "BHR", f"{SITE}_BHR")
    plot_QC(lastUpscale_dir, Date_albedo_upscale, BHRUNCupscale, Date_ALBEDO_gbov, BHRUNCgbov, "BHR~uncertainty", f"{SITE}_BHR_unc")
    plot_QC(lastUpscale_dir, Date_albedo_upscale, DHRupscale, Date_ALBEDO_gbov, DHRgbov, "DHR", f"{SITE}_DHR")
    plot_QC(lastUpscale_dir, Date_albedo_upscale, DHRUNCupscale, Date_ALBEDO_gbov, DHRUNCgbov, "DHR~uncertainty", f"{SITE}_DHR_unc")

    plot_output(lastUpscale_dir, Date_TOCR_upscale, TOCRupscale, TOCRUNCupscale, Date_TOCR_gbov, TOCRgbov, TOCRUNCgbov, "Top~of~canopy~reflectance~with~its~uncertainty", f"{SITE}_TOC-R_withUnc")
    plot_output(lastUpscale_dir, Date_albedo_upscale, BHRupscale, BHRUNCupscale, Date_ALBEDO_gbov, BHRgbov, BHRUNCgbov, "BHR~with~its~uncertainty", f"{SITE}_BHR_withUnc")
    plot_output(lastUpscale_dir, Date_albedo_upscale, DHRupscale, DHRUNCupscale, Date_ALBEDO_gbov, DHRgbov, DHRUNCgbov, "DHR~with~its~uncertainty", f"{SITE}_DHR_withUnc")

for i in range(len(Sites_list)):
    # Combine the coordinates and point labels
    data = {
        "Name": [f"Point {j + 1}" for j in range(9)],  # Point names Point 1, Point 2, ..., Point 9
        "Longitude": lon[i],
        "Latitude": lat[i]
    }

    # Create a DataFrame for the current site
    df = pd.DataFrame(data)

    # Save DataFrame to CSV file with site name
    filename = lastUpscale_dir + f"/QCheck/CSV/{Sites_list[i]}.csv"
    df.to_csv(filename, index=False)