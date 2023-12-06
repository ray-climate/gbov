#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_upscaling.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        05/12/2023 22:35

from datetime import datetime
from osgeo import gdal
import pandas as pd
import numpy as np
import os

# tower coordinate dic.
tower_coordinate = {}
tower_coordinate['Brasschaat'] = [51.3075, 4.5199, 'BRAS']

# instrument and canopy height.
canopy_height = {}
canopy_height['Brasschaat'] = [40., 21., 'BRAS']

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

def dhr_correction(sentinel2_dir, height_tower, height_canopy, lat, lon):

    radius = np.tan(np.deg2rad(85.)) * (height_tower - height_canopy)
    print('radius: ', radius)

    sentinel2_dhr_dir = os.path.join(sentinel2_dir, 'tile_0', 'albedo')
    for dhr_file in os.listdir(sentinel2_dhr_dir):
        if dhr_file.endswith('B02_UCL_bhr.jp2'):
            dhr_b02 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B03_UCL_bhr.jp2'):
            dhr_b03 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B04_UCL_bhr.jp2'):
            dhr_b04 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('BA8_UCL_bhr.jp2'):
            dhr_b8A = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B11_UCL_bhr.jp2'):
            dhr_b11 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))
        if dhr_file.endswith('B12_UCL_bhr.jp2'):
            dhr_b12 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))

    # find the dhr_b02 pixel within the radius, center at the given lat, lon
    geotransform = dhr_b02.GetGeoTransform()
    projection = dhr_b02.GetProjection()

    # Create a coordinate transformation object
    source_sr = gdal.osr.SpatialReference()
    source_sr.ImportFromWkt(projection)

    target_sr = gdal.osr.SpatialReference()
    target_sr.ImportFromEPSG(4326)  # EPSG code for WGS84

    coord_transform = gdal.osr.CoordinateTransformation(target_sr, source_sr)

    # Transform lat, lon to the raster's coordinate system
    tower_utm_x, tower_utm_y, _ = coord_transform.TransformPoint(lon, lat)

    # Raster size
    cols = dhr_b02.RasterXSize
    rows = dhr_b02.RasterYSize

    # Generate grid of pixel coordinates
    x_indices = np.arange(cols)
    y_indices = np.arange(rows)
    col_mesh, row_mesh = np.meshgrid(x_indices, y_indices)

    # Convert meshgrid to geospatial coordinates
    x_geo_mesh = geotransform[0] + col_mesh * geotransform[1] + row_mesh * geotransform[2]
    y_geo_mesh = geotransform[3] + col_mesh * geotransform[4] + row_mesh * geotransform[5]

    # Calculate the distance from the tower for all pixels
    distance_mesh = np.sqrt((x_geo_mesh - tower_utm_x) ** 2 + (y_geo_mesh - tower_utm_y) ** 2)
    print('distance_mesh: ', distance_mesh)
    # Find pixels within the specified radius
    pixels_within_radius = np.where(distance_mesh <= radius)
    print('pixels_within_radius: ', pixels_within_radius[0])
    quit()



def main():

    tower_retrieval_dir = '../ReferenceMeasurements/OUTPUT_dhr_bhr_tocr'
    sentinel2_dir = '/gws/nopw/j04/gbov/ruis/gbov/Sentinel2'

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

                dhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, lat, lon)
                quit()



if __name__ == '__main__':
    main()