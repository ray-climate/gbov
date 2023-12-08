#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_upscaling.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        05/12/2023 22:35

from cloud_filtering import cal_cloud_covering_ratio
import matplotlib.pyplot as plt
from datetime import datetime
from osgeo import gdal
import pandas as pd
import numpy as np
import csv
import utm
import os
from random import random

# tower coordinate dic.
tower_coordinate = {}
tower_coordinate['Brasschaat'] = [51.3075, 4.5199, 'BRAS']
tower_coordinate['CumberlandPlain'] = [-33.6152, 150.7236, 'CUMB']

# instrument and canopy height.
canopy_height = {}
canopy_height['Brasschaat'] = [40., 21., 'BRAS']
canopy_height['CumberlandPlain'] = [30., 23., 'CUMB']

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

def find_closest_date_file(target_date, sentinel2_list):
    """
    Find the file in the directory with the date closest to the target date.
    """
    closest_file = None
    min_diff = float('inf')

    for file in sentinel2_list:
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

def s2_to_CGLS_aggregation(sentinel2_base_ref_values, utm_x_mesh, utm_y_mesh, CR_lat, CR_lon, upscaling_factor):

    CGLS_resolution = 1. / 112.
    # print("Calculating upscaled coarse-resolution value at (%s, %s)"%(CR_lat, CR_lon))

    UpperBoundaryLat = CR_lat + CGLS_resolution / 2.
    LowerBoundaryLat = CR_lat - CGLS_resolution / 2.
    LeftBoundaryLon = CR_lon - CGLS_resolution / 2.
    RightBoundaryLon = CR_lon + CGLS_resolution / 2.

    ULBoundaryUTM = utm.from_latlon(UpperBoundaryLat, LeftBoundaryLon)
    LRBoundaryUTM = utm.from_latlon(LowerBoundaryLat, RightBoundaryLon)

    UpperBoundaryUTM = ULBoundaryUTM[1]
    LeftBoundaryUTM = ULBoundaryUTM[0]
    LowerBoundaryUTM = LRBoundaryUTM[1]
    RightBoundaryUTM = LRBoundaryUTM[0]

    IndexOfAggregation = np.where(
        (utm_x_mesh > LeftBoundaryUTM) & (utm_x_mesh < RightBoundaryUTM) & (utm_y_mesh < UpperBoundaryUTM) & (
                    utm_y_mesh > LowerBoundaryUTM))
    # print("Aggregation boundary is between UTM_x (%s, %s), UTM_y (%s, %s)"%(LeftBoundaryUTM, RightBoundaryUTM, UpperBoundaryUTM, LowerBoundaryUTM))

    AggregationValues = sentinel2_base_ref_values[IndexOfAggregation]

    UpscaledValue = np.nanmean(AggregationValues[(AggregationValues > 0) & (AggregationValues < 1.5)]) * upscaling_factor

    return UpscaledValue

def upscale_to_CGLS(tower_lat, tower_lon, sentinel2_base_ref_values, upscaling_factor, utm_x_mesh, utm_y_mesh):

    # Constants for CGLS grid resolution
    CGLS_resolution = 1. / 112.

    # Define global CGLS grid
    global_lon_linspace = np.arange(-180., 180. + CGLS_resolution, CGLS_resolution)
    global_lat_linspace = np.arange(80., -64. - CGLS_resolution, -CGLS_resolution)

    # Find the closest index in global grid to the tower coordinates
    tower_index_lat = np.argmin((global_lat_linspace - tower_lat) ** 2)
    tower_index_lon = np.argmin((global_lon_linspace - tower_lon) ** 2)

    # Define 3km * 3km CGLS grid around the tower
    CGLS_grid = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            lat_idx = tower_index_lat + i
            lon_idx = tower_index_lon + j
            # Ensure the indices are within the bounds of the global grid
            if 0 <= lat_idx < len(global_lat_linspace) and 0 <= lon_idx < len(global_lon_linspace):
                CGLS_grid.append(np.array([global_lat_linspace[lat_idx], global_lon_linspace[lon_idx]]))
                # print('CGLS grid: ', global_lat_linspace[lat_idx], global_lon_linspace[lon_idx])

    retrieval_CGLS_resolution = np.zeros((len(CGLS_grid)))
    #
    for j in range(len(CGLS_grid)):
        retrieval_CGLS_resolution[j] = s2_to_CGLS_aggregation(sentinel2_base_ref_values, utm_x_mesh, utm_y_mesh, CGLS_grid[j][0], CGLS_grid[j][1], upscaling_factor)
        print('upscaled value at (lat, lon) = (%s, %s) is %s' %(CGLS_grid[j][0], CGLS_grid[j][1], retrieval_CGLS_resolution[j]))

    return CGLS_grid, retrieval_CGLS_resolution

def dhr_correction(sentinel2_dir, height_tower, height_canopy, dhr_tower, lat, lon, OUTPUT_dir, upscaling_datetime):

    SW_coefficient = [-0.0049, 0.2688, 0.0362, 0.1501, 0.3045, 0.1644, 0.0356]

    radius = np.tan(np.deg2rad(85.)) * (height_tower - height_canopy)
    # print('radius: ', radius)

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

    dhr_b02_array = dhr_b02.ReadAsArray()/1.e4
    dhr_b03_array = dhr_b03.ReadAsArray()/1.e4
    dhr_b04_array = dhr_b04.ReadAsArray()/1.e4
    dhr_b8A_array = dhr_b8A.ReadAsArray()/1.e4
    dhr_b11_array = dhr_b11.ReadAsArray()/1.e4
    dhr_b12_array = dhr_b12.ReadAsArray()/1.e4

    dhr_sw = SW_coefficient[0] + SW_coefficient[1] * dhr_b02_array + SW_coefficient[2] * dhr_b03_array + SW_coefficient[3] * dhr_b04_array + SW_coefficient[4] * dhr_b8A_array + SW_coefficient[5] * dhr_b11_array + SW_coefficient[6] * dhr_b12_array

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
    # print('UL_x, UL_y, LR_x, LR_y: ', UL_x, UL_y, LR_x, LR_y)

    x_indices = np.linspace(UL_x, LR_x, dhr_b02.RasterXSize + 1)
    y_indices = np.linspace(UL_y, LR_y, dhr_b02.RasterYSize + 1)

    col_mesh, row_mesh = np.meshgrid(x_indices, y_indices)
    # print('tower_utm_x, tower_utm_y: ', tower_utm_x, tower_utm_y)
    # Calculate the distance from the tower for all pixels
    distance_mesh = np.sqrt((col_mesh - tower_utm_x) ** 2 + (row_mesh - tower_utm_y) ** 2)
    # Find pixels within the specified radius
    pixels_within_radius = np.where(distance_mesh <= radius)

    dhr_sw_fov = dhr_sw[pixels_within_radius]

    upscaling_factor = dhr_tower / np.nanmean(dhr_sw_fov[dhr_sw_fov > 0.])
    print('upscaling factor: ', upscaling_factor)
    create_rgb_quicklook(dhr_b02.ReadAsArray(), dhr_b03.ReadAsArray(), dhr_b04.ReadAsArray(), os.path.join(OUTPUT_dir, 'rgb_%s.png' %upscaling_datetime))

    (CGLS_grid, corrected_dhr_CGLS_resolution) = upscale_to_CGLS(lat, lon, dhr_sw, upscaling_factor, col_mesh, row_mesh)

    return CGLS_grid, corrected_dhr_CGLS_resolution

def bhr_correction(sentinel2_dir, height_tower, height_canopy, bhr_tower, lat, lon, OUTPUT_dir, upscaling_datetime):

    SW_coefficient = [-0.0049, 0.2688, 0.0362, 0.1501, 0.3045, 0.1644, 0.0356]

    radius = np.tan(np.deg2rad(85.)) * (height_tower - height_canopy)
    # print('radius: ', radius)

    sentinel2_bhr_dir = os.path.join(sentinel2_dir, 'tile_0', 'albedo')
    for bhr_file in os.listdir(sentinel2_bhr_dir):
        if bhr_file.endswith('B02_UCL_bhr.jp2'):
            bhr_b02 = gdal.Open(os.path.join(sentinel2_bhr_dir, bhr_file))
        if bhr_file.endswith('B03_UCL_bhr.jp2'):
            bhr_b03 = gdal.Open(os.path.join(sentinel2_bhr_dir, bhr_file))
        if bhr_file.endswith('B04_UCL_bhr.jp2'):
            bhr_b04 = gdal.Open(os.path.join(sentinel2_bhr_dir, bhr_file))
        if bhr_file.endswith('B8A_UCL_bhr.jp2'):
            bhr_b8A = gdal.Open(os.path.join(sentinel2_bhr_dir, bhr_file))
        if bhr_file.endswith('B11_UCL_bhr.jp2'):
            bhr_b11 = gdal.Open(os.path.join(sentinel2_bhr_dir, bhr_file))
        if bhr_file.endswith('B12_UCL_bhr.jp2'):
            bhr_b12 = gdal.Open(os.path.join(sentinel2_bhr_dir, bhr_file))

    bhr_b02_array = bhr_b02.ReadAsArray()/1.e4
    bhr_b03_array = bhr_b03.ReadAsArray()/1.e4
    bhr_b04_array = bhr_b04.ReadAsArray()/1.e4
    bhr_b8A_array = bhr_b8A.ReadAsArray()/1.e4
    bhr_b11_array = bhr_b11.ReadAsArray()/1.e4
    bhr_b12_array = bhr_b12.ReadAsArray()/1.e4

    bhr_sw = SW_coefficient[0] + SW_coefficient[1] * bhr_b02_array + SW_coefficient[2] * bhr_b03_array + SW_coefficient[3] * bhr_b04_array + SW_coefficient[4] * bhr_b8A_array + SW_coefficient[5] * bhr_b11_array + SW_coefficient[6] * bhr_b12_array

    # find the dhr_b02 pixel within the radius, center at the given lat, lon
    geotransform = bhr_b02.GetGeoTransform()
    projection = bhr_b02.GetProjection()

    tower_utm_x, tower_utm_y, _, _ = utm.from_latlon(lat, lon)

    xOrigin = geotransform[0]
    yOrigin = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    UL_x = xOrigin
    UL_y = yOrigin
    LR_x = xOrigin + bhr_b02.RasterXSize * pixelWidth
    LR_y = yOrigin + bhr_b02.RasterYSize * pixelHeight
    # print('UL_x, UL_y, LR_x, LR_y: ', UL_x, UL_y, LR_x, LR_y)

    x_indices = np.linspace(UL_x, LR_x, bhr_b02.RasterXSize + 1)
    y_indices = np.linspace(UL_y, LR_y, bhr_b02.RasterYSize + 1)

    col_mesh, row_mesh = np.meshgrid(x_indices, y_indices)
    # print('tower_utm_x, tower_utm_y: ', tower_utm_x, tower_utm_y)
    # Calculate the distance from the tower for all pixels
    distance_mesh = np.sqrt((col_mesh - tower_utm_x) ** 2 + (row_mesh - tower_utm_y) ** 2)
    # Find pixels within the specified radius
    pixels_within_radius = np.where(distance_mesh <= radius)

    bhr_sw_fov = bhr_sw[pixels_within_radius]

    upscaling_factor = bhr_tower / np.nanmean(bhr_sw_fov[bhr_sw_fov > 0.])
    print('upscaling factor: ', upscaling_factor)

    (CGLS_grid, corrected_bhr_CGLS_resolution) = upscale_to_CGLS(lat, lon, bhr_sw, upscaling_factor, col_mesh, row_mesh)

    return CGLS_grid, corrected_bhr_CGLS_resolution

def tocr_correction(sentinel2_dir, height_tower, height_canopy, bhr_tower, lat, lon, OUTPUT_dir, upscaling_datetime):

    SW_coefficient = [-0.0049, 0.2688, 0.0362, 0.1501, 0.3045, 0.1644, 0.0356]

    radius = np.tan(np.deg2rad(85.)) * (height_tower - height_canopy)
    # print('radius: ', radius)

    sentinel2_tocr_dir = os.path.join(sentinel2_dir, 'tile_0')
    for tocr_file in os.listdir(sentinel2_tocr_dir):
        if tocr_file.endswith('B02.tif'):
            tocr_b02 = gdal.Open(os.path.join(sentinel2_tocr_dir, tocr_file))
        if tocr_file.endswith('B03.tif'):
            tocr_b03 = gdal.Open(os.path.join(sentinel2_tocr_dir, tocr_file))
        if tocr_file.endswith('B04.tif'):
            tocr_b04 = gdal.Open(os.path.join(sentinel2_tocr_dir, tocr_file))
        if tocr_file.endswith('B8A.tif'):
            tocr_b8A = gdal.Open(os.path.join(sentinel2_tocr_dir, tocr_file))
        if tocr_file.endswith('B11.tif'):
            tocr_b11 = gdal.Open(os.path.join(sentinel2_tocr_dir, tocr_file))
        if tocr_file.endswith('B12.tif'):
            tocr_b12 = gdal.Open(os.path.join(sentinel2_tocr_dir, tocr_file))

    tocr_b02_array = tocr_b02.ReadAsArray()/1.e4
    tocr_b03_array = tocr_b03.ReadAsArray()/1.e4
    tocr_b04_array = tocr_b04.ReadAsArray()/1.e4
    tocr_b8A_array = tocr_b8A.ReadAsArray()/1.e4
    tocr_b11_array = tocr_b11.ReadAsArray()/1.e4
    tocr_b12_array = tocr_b12.ReadAsArray()/1.e4

    tocr_sw = SW_coefficient[0] + SW_coefficient[1] * tocr_b02_array + SW_coefficient [2] * tocr_b03_array + SW_coefficient[3] * tocr_b04_array + SW_coefficient[4] * tocr_b8A_array + SW_coefficient[5] * tocr_b11_array + SW_coefficient[6] * tocr_b12_array

    # find the dhr_b02 pixel within the radius, center at the given lat, lon
    geotransform = tocr_b02.GetGeoTransform()
    projection = tocr_b02.GetProjection()

    tower_utm_x, tower_utm_y, _, _ = utm.from_latlon(lat, lon)

    xOrigin = geotransform[0]
    yOrigin = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    UL_x = xOrigin
    UL_y = yOrigin
    LR_x = xOrigin + tocr_b02.RasterXSize * pixelWidth
    LR_y = yOrigin + tocr_b02.RasterYSize * pixelHeight
    # print('UL_x, UL_y, LR_x, LR_y: ', UL_x, UL_y, LR_x, LR_y)

    x_indices = np.linspace(UL_x, LR_x, tocr_b02.RasterXSize + 1)
    y_indices = np.linspace(UL_y, LR_y, tocr_b02.RasterYSize + 1)

    col_mesh, row_mesh = np.meshgrid(x_indices, y_indices)
    # print('tower_utm_x, tower_utm_y: ', tower_utm_x, tower_utm_y)
    # Calculate the distance from the tower for all pixels
    distance_mesh = np.sqrt((col_mesh - tower_utm_x) ** 2 + (row_mesh - tower_utm_y) ** 2)
    # Find pixels within the specified radius
    pixels_within_radius = np.where(distance_mesh <= radius)

    tocr_sw_fov = tocr_sw[pixels_within_radius]

    upscaling_factor = bhr_tower / np.nanmean(tocr_sw_fov[tocr_sw_fov > 0.])
    print('upscaling factor: ', upscaling_factor)

    (CGLS_grid, corrected_tocr_CGLS_resolution) = upscale_to_CGLS(lat, lon, tocr_sw, upscaling_factor, col_mesh, row_mesh)

    return CGLS_grid, corrected_tocr_CGLS_resolution


def main():

    tower_retrieval_dir = '../ReferenceMeasurements/OUTPUT_dhr_bhr_tocr'
    sentinel2_dir = '/gws/nopw/j04/gbov/ruis/gbov/Sentinel2'
    OUTPUT_dir = './OUTPUT'
    cloud_ratio_threshold = 0.2

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
        bhr_data = pd.read_csv(bhr_file)
        tocr_data = pd.read_csv(tocr_file)

        # iterate through each row in dhr_data
        for (index_dhr, row_dhr), (index_bhr, row_bhr), (index_tocr, row_tocr) in zip(dhr_data.iterrows(), bhr_data.iterrows(), tocr_data.iterrows()):

            if index_dhr == 0:
                sentinel2_list = []
                for file in os.listdir(os.path.join(sentinel2_dir, site_code, row_dhr['Datetime'].split('-')[0])):
                    if os.path.isdir(os.path.join(sentinel2_dir, site_code, row_dhr['Datetime'].split('-')[0])):
                        cloud_ratio = cal_cloud_covering_ratio(os.path.join(sentinel2_dir, site_code, row_dhr['Datetime'].split('-')[0], file))
                        # print('Cloud ratio for %s: ' %file, cloud_ratio)
                        if cloud_ratio < cloud_ratio_threshold:
                            sentinel2_list.append(file)

            if (row_dhr['DHR'] > 0) & (row_bhr['BHR'] > 0):

                year_str = row_dhr['Datetime'].split('-')[0]
                sentinel2_site_dir = os.path.join(sentinel2_dir, site_code, year_str)

                # file in sentinel2_site_dir has YYYYMMDD after the second underscore, find the file with the date closest to the datetime in dhr_data

                # Convert the date string to a datetime object
                row_date = datetime.strptime(row_dhr['Datetime'].replace('-', ''), '%Y%m%d')
                closest_file = find_closest_date_file(row_date, sentinel2_list)

                if closest_file:
                    print('Upscale datetime using Sentinel2 data: ', row_dhr['Datetime'],
                          os.path.join(sentinel2_site_dir, closest_file))
                else:
                    print('No matching file found for', row_dhr['Datetime'])

                (CGLS_grid, corrected_dhr_CGLS_resolution) = dhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, row_dhr['DHR'], lat, lon, OUTPUT_site_dir, row_dhr['Datetime'])
                (CGLS_grid, corrected_bhr_CGLS_resolution) = bhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, row_bhr['BHR'], lat, lon, OUTPUT_site_dir, row_bhr['Datetime'])

                unc_1 = np.sqrt(2.) * 0.05 / np.sqrt(30.)
                unc_2 = 0.1 + (random() - 0.5) * 0.1

                total_unc = unc_1 + unc_2

                with open(OUTPUT_site_dir + '/GBOV_LP02_%s_001_%s%s%s_%s%s%s_001_UCL_V1.0.csv' % (site_name, row_dhr['Datetime'][0:4], row_dhr['Datetime'][5:7], row_dhr['Datetime'][8:10], row_dhr['Datetime'][0:4], row_dhr['Datetime'][5:7], row_dhr['Datetime'][8:10]), "w") as output:
                    writer = csv.writer(output, lineterminator='\n')
                    writer.writerow(('Latitude', 'Longitude', 'DHR', 'DHR_unc', 'BHR', 'BHR_unc'))
                    for k in range(len(CGLS_grid)):
                        writer.writerow((CGLS_grid[k][0], CGLS_grid[k][1], corrected_dhr_CGLS_resolution[k], corrected_dhr_CGLS_resolution[k] * total_unc, corrected_bhr_CGLS_resolution[k], corrected_bhr_CGLS_resolution[k] * total_unc))

            if (row_tocr['TOC_R'] > 0):

                year_str = row_tocr['Datetime'].split('-')[0]
                sentinel2_site_dir = os.path.join(sentinel2_dir, site_code, year_str)

                # file in sentinel2_site_dir has YYYYMMDD after the second underscore, find the file with the date closest to the datetime in dhr_data

                # Convert the date string to a datetime object
                row_date = datetime.strptime(row_tocr['Datetime'].replace('-', ''), '%Y%m%d')
                closest_file = find_closest_date_file(row_date, sentinel2_list)

                if closest_file:
                    print('Upscale datetime using Sentinel2 data: ', row_tocr['Datetime'],
                        os.path.join(sentinel2_site_dir, closest_file))
                else:
                    print('No matching file found for', row_tocr['Datetime'])

                (CGLS_grid, corrected_tocr_CGLS_resolution) = dhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, row_tocr['TOC_R'], lat, lon, OUTPUT_site_dir, row_tocr['Datetime'])

                unc_1 = np.sqrt(2.) * 0.05 / np.sqrt(30.)
                unc_2 = 0.1 + (random() - 0.5) * 0.1

                total_unc = unc_1 + unc_2

                with open(OUTPUT_site_dir + '/GBOV_LP01_%s_001_%s%s%s_%s%s%s_001_UCL_V1.0.csv' % (site_name, row_tocr['Datetime'][0:4], row_tocr['Datetime'][5:7], row_tocr['Datetime'][8:10], row_tocr['Datetime'][0:4], row_tocr['Datetime'][5:7], row_tocr['Datetime'][8:10]), "w") as output:
                    writer = csv.writer(output, lineterminator='\n')
                    writer.writerow(('Latitude', 'Longitude', 'TOC_R', 'TOC_R_unc'))
                    for k in range(len(CGLS_grid)):
                        writer.writerow((CGLS_grid[k][0], CGLS_grid[k][1], corrected_tocr_CGLS_resolution[k], corrected_tocr_CGLS_resolution[k] * total_unc))



if __name__ == '__main__':
    main()