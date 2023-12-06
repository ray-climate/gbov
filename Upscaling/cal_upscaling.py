#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_upscaling.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        05/12/2023 22:35

import matplotlib.pyplot as plt
from datetime import datetime
from osgeo import gdal
import pandas as pd
import numpy as np
import utm
import os

# tower coordinate dic.
tower_coordinate = {}
tower_coordinate['Brasschaat'] = [51.3075, 4.5199, 'BRAS']

# instrument and canopy height.
canopy_height = {}
canopy_height['Brasschaat'] = [40., 21., 'BRAS']

# a function to generate RGB quicklook image for each site using blue (B2), green (B3), red (B4) bands.
def create_rgb_quicklook(band2, band3, band4, output_file):
    def normalize(array):
        # Replace negative values with NaN to handle separately
        array = np.where(array < 0, np.nan, array)
        array_min, array_max = np.nanmin(array), np.nanmax(array)
        # Normalize the array, ignoring NaNs
        return (array - array_min) / (array_max - array_min)

    # Normalize the bands
    band_red = normalize(band4)
    band_green = normalize(band3)
    band_blue = normalize(band2)

    # Create an RGB array
    rgb = np.dstack((band_red, band_green, band_blue))

    # Replace NaNs (originally negative values) with grey
    nan_mask = np.isnan(rgb)
    rgb[nan_mask] = 0.5  # Set grey value (0.5 in normalized scale)

    # Adjust brightness
    brightness_factor = 1.5  # Adjust this factor to control brightness
    rgb = np.clip(rgb * brightness_factor, 0, 1)

    # Save the image
    plt.imsave(output_file, rgb, cmap='gray')


def find_closest_date_file(target_date, directory):
    """
    Find the file in the directory with the date closest to the target date.
    """
    closest_file = None
    min_diff = float('inf')

    for file in os.listdir(directory):
        # Extracting the date from the file name (assuming YYYYMMDD format after second underscore)
        try:
            file_date_str = file.split('_')[2]
            file_date = datetime.strptime(file_date_str, '%Y%m%d')
            diff = abs((target_date - file_date).days)

            if diff < min_diff:
                min_diff = diff
                closest_file = file
        except (IndexError, ValueError):
            # Handle files that do not match the expected format
            continue

    return closest_file

def dhr_correction(sentinel2_dir, height_tower, height_canopy, lat, lon, OUTPUT_dir, upscaling_datetime):

    radius = np.tan(np.deg2rad(85.)) * (height_tower - height_canopy)
    print('radius: ', radius)

    sentinel2_dhr_dir = os.path.join(sentinel2_dir, 'tile_0', 'albedo')
    for dhr_file in os.listdir(sentinel2_dhr_dir):
        if dhr_file.endswith('B02_UCL_dhr.jp2'):
            dhr_b02 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B03_UCL_dhr.jp2'):
            dhr_b03 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B04_UCL_dhr.jp2'):
            dhr_b04 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B8A_UCL_dhr.jp2'):
            dhr_b8A = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B11_UCL_dhr.jp2'):
            dhr_b11 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B12_UCL_dhr.jp2'):
            dhr_b12 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))

    # find the dhr_b02 pixel within the radius, center at the given lat, lon
    geotransform = dhr_b02.GetGeoTransform()
    projection = dhr_b02.GetProjection()

    tower_utm_x, tower_utm_y, _, _ = utm.from_latlon(lat, lon)

    xOrigin = geotransform[0]
    yOrigin = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    UL_x = xOrigin
    UL_y = yOrigin
    LR_x = xOrigin + dhr_b02.RasterXSize * pixelWidth
    LR_y = yOrigin + dhr_b02.RasterYSize * pixelHeight
    print('UL_x, UL_y, LR_x, LR_y: ', UL_x, UL_y, LR_x, LR_y)

    x_indices = np.linspace(UL_x, LR_x, dhr_b02.RasterXSize + 1)
    y_indices = np.linspace(UL_y, LR_y, dhr_b02.RasterYSize + 1)

    col_mesh, row_mesh = np.meshgrid(x_indices, y_indices)
    print('tower_utm_x, tower_utm_y: ', tower_utm_x, tower_utm_y)
    # Calculate the distance from the tower for all pixels
    distance_mesh = np.sqrt((col_mesh - tower_utm_x) ** 2 + (row_mesh - tower_utm_y) ** 2)
    print('distance_mesh: ', distance_mesh)
    # Find pixels within the specified radius
    pixels_within_radius = np.where(distance_mesh <= radius)
    # get the last part of filename from sentinel2_dhr_dir
    filename = sentinel2_dhr_dir.split('/')[-3]
    print('filename: ', filename)
    quit()

    create_rgb_quicklook(dhr_b02.ReadAsArray(), dhr_b03.ReadAsArray(), dhr_b04.ReadAsArray(), os.path.join(OUTPUT_dir, 'rgb_%s.png' %upscaling_datetime))
    quit()





def main():

    tower_retrieval_dir = '../ReferenceMeasurements/OUTPUT_dhr_bhr_tocr'
    sentinel2_dir = '/gws/nopw/j04/gbov/ruis/gbov/Sentinel2'
    OUTPUT_dir = './OUTPUT'

    os.makedirs(OUTPUT_dir, exist_ok=True)

    upscale_filelist = []
    for file in os.listdir(tower_retrieval_dir):
        if file.endswith('GBOV_RM01_Brasschaat_001_20210101T000000Z_20211231T233000Z_007_ACR_V2.0_DHR.csv'):
            upscale_filelist.append(file[:-8])
            print('Site to be upscaled: ', file[:-8])

    for tower_file in upscale_filelist:

        dhr_file = os.path.join(tower_retrieval_dir, tower_file + '_DHR.csv')
        bhr_file = os.path.join(tower_retrieval_dir, tower_file + '_BHR.csv')
        tocr_file = os.path.join(tower_retrieval_dir, tower_file + '_TOC_R.csv')

        # extract site name from file name
        site_name = tower_file.split('_')[2]
        site_code = tower_coordinate[site_name][2]
        print('Site name, site code: ', site_name, site_code)

        OUTPUT_site_dir = os.path.join(OUTPUT_dir, site_name)
        os.makedirs(OUTPUT_site_dir, exist_ok=True)

        lat, lon = tower_coordinate[site_name][0], tower_coordinate[site_name][1]
        height_tower, height_canopy = canopy_height[site_name][0], canopy_height[site_name][1]
        print('Tower height, canopy height: ', height_tower, height_canopy)

        # read values in dhr_file
        dhr_data = pd.read_csv(dhr_file)
        # iterate through each row in dhr_data
        for index, row in dhr_data.iterrows():
            # print(row['datetime'], row['dhr'])
            if row['DHR'] > 0:
                year_str = row['Datetime'].split('-')[0]
                sentinel2_site_dir = os.path.join(sentinel2_dir, site_code, year_str)

                # file in sentinel2_site_dir has YYYYMMDD after the second underscore, find the file with the date closest to the datetime in dhr_data

                # Convert the date string to a datetime object
                row_date = datetime.strptime(row['Datetime'].replace('-', ''), '%Y%m%d')
                closest_file = find_closest_date_file(row_date, sentinel2_site_dir)

                if closest_file:
                    print('Upscale datetime using Sentinel2 data: ', row['Datetime'],
                          os.path.join(sentinel2_site_dir, closest_file))
                else:
                    print('No matching file found for', row['Datetime'])

                dhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, lat, lon, OUTPUT_site_dir, row['Datetime'])
                quit()


if __name__ == '__main__':
    main()