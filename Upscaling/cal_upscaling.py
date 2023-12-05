#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_upscaling.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        05/12/2023 22:35

import pandas as pd
import os

def main():

    rm1_dir = '../ReferenceMeasurements/RM1_data'
    tower_retrieval_dir = '../ReferenceMeasurements/OUTPUT_dhr_bhr_tocr'

    upscale_filelist = []
    for file in os.listdir(rm1_dir):
        if file.endswith('GBOV_RM01_Brasschaat_001_20210101T000000Z_20211231T233000Z_007_ACR_V2.0.csv'):
            upscale_filelist.append(file[:-4])
            print('Site to be upscaled: ', file[:-4])

    for file in upscale_filelist:

        dhr_file = os.path.join(tower_retrieval_dir, file + '_DHR.csv')
        bhr_file = os.path.join(tower_retrieval_dir, file + '_BHR.csv')
        tocr_file = os.path.join(tower_retrieval_dir, file + '_TOC_R.csv')

        # read values in dhr_file
        dhr_data = pd.read_csv(dhr_file)
        # iterate through each row in dhr_data
        for index, row in dhr_data.iterrows():
            # print(row['datetime'], row['dhr'])
            if row['DHR'] > 0:
                print(row['Datetime'], row['DHR'])

if __name__ == '__main__':
    main()