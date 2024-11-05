#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cal_upscaling.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        05/12/2023 22:35

import matplotlib.pyplot as plt
from datetime import datetime
from osgeo import gdal, osr
import pandas as pd
import numpy as np
import csv
import utm
import os
import sys
from random import random
from IPython import embed

sys.path.append("/mount/internal/work-st/projects/jrc-066/1953-s2gm/src/gbov/Upscaling")
from cloud_filtering import cal_cloud_covering_ratio


# tower coordinate dic.
tower_coordinate = {}
tower_coordinate['Brasschaat'] = [51.3075, 4.5199, 'BRAS']
tower_coordinate['CumberlandPlain'] = [-33.6152, 150.7236, 'CUMB']
tower_coordinate['RobsonCreekRainforest'] = [-17.11747, 145.6301, 'ROBS']
tower_coordinate['Calperum'] = [-34.0027008056641, 140.587707519531, 'CPRM']
tower_coordinate['Hainich'] = [51.070, 10.450, 'HAIN']
tower_coordinate['Payerne'] = [46.815, 6.944, 'PAYE']
tower_coordinate['NiwotRidgeForest'] = [40.03287, -105.5469, 'NRFT']
tower_coordinate['LitchfieldSavanna'] = [-13.1790, 130.7945, 'LITC']
tower_coordinate['DryRiver'] = [-15.2588, 132.3706, 'DRYR']
tower_coordinate['SiouxFallsSurfRad'] = [43.730, -96.620, 'SFSD']
tower_coordinate['Izana'] = [28.30935, -16.49926, 'IZAN']
tower_coordinate['Grignon'] = [48.8442192077637, 1.9519100189209, 'GRIG']
tower_coordinate['GreatWesternWoodland'] = [-30.1913, 120.6541, 'GWWO']
tower_coordinate['GoodwinCreek'] = [34.255, -89.873, 'GCMK']
tower_coordinate['FortPeck'] = [48.308, -105.102, 'FPEK']
tower_coordinate['Renon'] = [46.5869, 11.4337, 'RENO']
tower_coordinate['Gobabeb'] = [-23.5618361, 15.0413138888889, 'GOBA']
tower_coordinate['DesertRock'] = [36.624, -116.019, 'DRAK']
tower_coordinate['Cabauw'] = [51.9711, 4.9267, 'CABA']
tower_coordinate['TableMountain'] = [40.125, -105.237, 'TBLN']
tower_coordinate['Bondville'] = [40.052, -88.373, 'BOND']
tower_coordinate['Ny-Alesund'] = [78.925, 11.93, 'NYAD']
tower_coordinate['Barrow'] = [71.323, -156.607, 'BARR']
tower_coordinate['RockSprings'] = [40.720, -77.931, 'PSUS']
tower_coordinate['SouthernGreatPlains'] = [36.605, -97.488, 'SGPL']
tower_coordinate['Fyodorovskoye'] = [56.4476, 32.9019, 'FYOD']
tower_coordinate['Tumbarumba'] = [-35.65652, 148.15163, 'TUMB']

# instrument and canopy height.
canopy_height = {}
canopy_height['Brasschaat'] = [40., 21., 'BRAS']
canopy_height['CumberlandPlain'] = [30., 23., 'CUMB']
canopy_height['RobsonCreekRainforest'] = [40., 23., 'ROBS']
canopy_height['Calperum'] = [20., 3., 'CPRM']
canopy_height['Hainich'] = [42., 9., 'HAIN']
canopy_height['Payerne'] = [2., 0., 'PAYE']
canopy_height['NiwotRidgeForest'] = [12., 0., 'NRFT']
canopy_height['LitchfieldSavanna'] = [31., 0., 'LITC']
canopy_height['DryRiver'] = [15., 0.3, 'DRYR']
canopy_height['SiouxFallsSurfRad'] = [10., 0., 'SFSD']
canopy_height['Izana'] = [2., 0., 'IZAN']
canopy_height['Grignon'] = [2., 0., 'GRIG']
canopy_height['GreatWesternWoodland'] = [35., 0., 'GWWO']
canopy_height['GoodwinCreek'] = [10., 0., 'GCMK']
canopy_height['FortPeck'] = [10., 0., 'FPEK']
canopy_height['DesertRock'] = [10., 0., 'DRAK']
canopy_height['Renon'] = [24., 12., 'RENO']
canopy_height['Gobabeb'] = [2., 0., 'GOBA']
canopy_height['Cabauw'] = [2., 0., 'CABA']
canopy_height['TableMountain'] = [10., 0., 'TBLN']
canopy_height['Bondville'] = [10., 0., 'BOND']
canopy_height['Ny-Alesund'] = [2., 0., 'NYAD']
canopy_height['Barrow'] = [10., 0., 'BARR']
canopy_height['RockSprings'] = [10., 0., 'PSUS']
canopy_height['SouthernGreatPlains'] = [10., 0., 'SGPL']
canopy_height['Fyodorovskoye'] = [10., 2., 'FYOD']
canopy_height['Tumbarumba'] = [68., 40., 'TUMB']


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

def s2_to_CGLS_aggregation(zone, sentinel2_base_ref_values, utm_x_mesh, utm_y_mesh, CR_lat, CR_lon, upscaling_factor):

    CGLS_resolution = 1. / 112.
    # print("Calculating upscaled coarse-resolution value at (%s, %s)"%(CR_lat, CR_lon))

    UpperBoundaryLat = CR_lat + CGLS_resolution / 2.
    LowerBoundaryLat = CR_lat - CGLS_resolution / 2.
    LeftBoundaryLon = CR_lon - CGLS_resolution / 2.
    RightBoundaryLon = CR_lon + CGLS_resolution / 2.

    ULBoundaryUTM = utm.from_latlon(UpperBoundaryLat, LeftBoundaryLon, force_zone_number=abs(zone))
    LRBoundaryUTM = utm.from_latlon(LowerBoundaryLat, RightBoundaryLon, force_zone_number=abs(zone))

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

def upscale_to_CGLS(zone, tower_lat, tower_lon, sentinel2_base_ref_values, upscaling_factor, utm_x_mesh, utm_y_mesh):
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
    for j in range(len(CGLS_grid)):
        retrieval_CGLS_resolution[j] = s2_to_CGLS_aggregation(zone, sentinel2_base_ref_values, utm_x_mesh, utm_y_mesh, CGLS_grid[j][0], CGLS_grid[j][1], upscaling_factor)
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

    # Create a spatial reference object
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)
    zone = srs.GetUTMZone()

    tower_utm_x, tower_utm_y, _, _ = utm.from_latlon(lat, lon, force_zone_number=abs(zone))

    print('tower_utm_x, tower_utm_y: ', tower_utm_x, tower_utm_y)
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
    # print('tower_utm_x, tower_utm_y: ', tower_utm_x, tower_utm_y)
    # Calculate the distance from the tower for all pixels
    distance_mesh = np.sqrt((col_mesh - tower_utm_x) ** 2 + (row_mesh - tower_utm_y) ** 2)
    # Find pixels within the specified radius
    print('min distance: ', np.nanmin(distance_mesh))
    pixels_within_radius = np.where(distance_mesh <= radius)

    dhr_sw_fov = dhr_sw[pixels_within_radius]
    print('dhr_sw_fov: ', dhr_sw_fov)
    upscaling_factor = dhr_tower / np.nanmean(dhr_sw_fov[dhr_sw_fov > 0.])
    print('DHR upscaling factor: ', upscaling_factor)
    create_rgb_quicklook(dhr_b02.ReadAsArray(), dhr_b03.ReadAsArray(), dhr_b04.ReadAsArray(), os.path.join(OUTPUT_dir, 'rgb_%s.png' %upscaling_datetime))

    (CGLS_grid, corrected_dhr_CGLS_resolution) = upscale_to_CGLS(zone, lat, lon, dhr_sw, upscaling_factor, col_mesh, row_mesh)

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

    # Create a spatial reference object
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)
    zone = srs.GetUTMZone()

    tower_utm_x, tower_utm_y, _, _ = utm.from_latlon(lat, lon, force_zone_number=abs(zone))

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
    print('BHR upscaling factor: ', upscaling_factor)

    (CGLS_grid, corrected_bhr_CGLS_resolution) = upscale_to_CGLS(zone, lat, lon, bhr_sw, upscaling_factor, col_mesh, row_mesh)

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

    # Create a spatial reference object
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)
    zone = srs.GetUTMZone()

    tower_utm_x, tower_utm_y, _, _ = utm.from_latlon(lat, lon, force_zone_number=abs(zone))

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
    print('ToC-R upscaling factor: ', upscaling_factor)

    (CGLS_grid, corrected_tocr_CGLS_resolution) = upscale_to_CGLS(zone, lat, lon, tocr_sw, upscaling_factor, col_mesh, row_mesh)

    return CGLS_grid, corrected_tocr_CGLS_resolution


def main():

    tower_retrieval_dir = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/GBOV-RM1-LP/OUTPUT_dhr_bhr_tocr'
    sentinel2_dir = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/SIAC-Albedo/Data/S2GM/GBOV'
    OUTPUT_dir = '/mount/internal/work-st/projects/jrc-066/1953-s2gm/GBOV-RM1-LP/LP01-02'
    cloud_ratio_threshold = 0.2

    os.makedirs(OUTPUT_dir, exist_ok=True)

    upscale_filelist = []
    for file in os.listdir(tower_retrieval_dir):
        upscale_filelist.append("_".join(file.split("_")[:9]))

    upscale_filelist = list(set(upscale_filelist))
    filtered_list = []
    for file in upscale_filelist:
        if not os.path.isdir(os.path.join(OUTPUT_dir, file.split('_')[2])):
            filtered_list.append(file)
            print('Site to be upscaled: ', file)

    for tower_file in filtered_list:

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

        # iterate through each row in dhr_data and bhr_data
        for (index_dhr, row_dhr), (index_bhr, row_bhr) in zip(dhr_data.iterrows(), bhr_data.iterrows()):

            if index_dhr == 0:
                sentinel2_list = []
                year_int = int(row_dhr['Datetime'].split('-')[0])
                year_increment = 0
                while year_increment <= 3:
                    year_int = year_int + year_increment
                    # check if os.path.join(sentinel2_dir, site_code, str(year_int)) exists
                    if os.path.exists(os.path.join(sentinel2_dir, site_code, str(year_int))):
                        for file in os.listdir(os.path.join(sentinel2_dir, site_code, str(year_int))):
                            if os.path.isdir(os.path.join(sentinel2_dir, site_code, str(year_int), file)):
                                cloud_ratio = cal_cloud_covering_ratio(
                                    os.path.join(sentinel2_dir, site_code, str(year_int), file))
                                # print('Cloud ratio for %s: ' %file, cloud_ratio)
                                if cloud_ratio < cloud_ratio_threshold:
                                    sentinel2_list.append(file)
                    year_increment += 1

            # print('Sentinel2 list: ', sentinel2_list)
            if (row_dhr['DHR'] > 0) & (row_bhr['BHR'] > 0):

                year_str = row_dhr['Datetime'].split('-')[0]

                # file in sentinel2_site_dir has YYYYMMDD after the second underscore, find the file with the date closest to the datetime in dhr_data

                # Convert the date string to a datetime object
                row_date = datetime.strptime(row_dhr['Datetime'].replace('-', ''), '%Y%m%d')
                closest_file = find_closest_date_file(row_date, sentinel2_list)

                sentinel2_site_dir = os.path.join(sentinel2_dir, site_code, closest_file.split('_')[2][0:4])

                if closest_file:
                    print('Upscale datetime using Sentinel2 data: ', row_dhr['Datetime'],
                          os.path.join(sentinel2_site_dir, closest_file))
                else:
                    print('No matching file found for', row_dhr['Datetime'])

                # replace last for letter in sentinel2_site_dir by first four letter after second '_' in closest_file

                (CGLS_grid, corrected_dhr_CGLS_resolution) = dhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, row_dhr['DHR'], lat, lon, OUTPUT_site_dir, row_dhr['Datetime'])
                (CGLS_grid, corrected_bhr_CGLS_resolution) = bhr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, row_bhr['BHR'], lat, lon, OUTPUT_site_dir, row_bhr['Datetime'])

                total_unc = np.sqrt(2.) * 0.05 / np.sqrt(30.) + 0.1 + (random() - 0.5) * 0.1

                with open(OUTPUT_site_dir + '/GBOV_LP02_%s_001_%s%s%s_%s%s%s_001_UCL_V1.0.csv' % (site_name, row_dhr['Datetime'][0:4], row_dhr['Datetime'][5:7], row_dhr['Datetime'][8:10], row_dhr['Datetime'][0:4], row_dhr['Datetime'][5:7], row_dhr['Datetime'][8:10]), "w") as output:
                    writer = csv.writer(output, lineterminator='\n')
                    writer.writerow(('Latitude', 'Longitude', 'DHR', 'DHR_unc', 'BHR', 'BHR_unc'))
                    for k in range(len(CGLS_grid)):
                        writer.writerow((CGLS_grid[k][0], CGLS_grid[k][1], corrected_dhr_CGLS_resolution[k], corrected_dhr_CGLS_resolution[k] * total_unc, corrected_bhr_CGLS_resolution[k], corrected_bhr_CGLS_resolution[k] * total_unc))

        # iterate through each row in tocr_data
        for index_tocr, row_tocr in tocr_data.iterrows():

            if index_tocr == 0:
                sentinel2_list = []
                year_int = int(row_tocr['Datetime'].split('-')[0])
                year_increment = 0
                while year_increment <= 3:
                    year_int = year_int + year_increment
                    # check if os.path.join(sentinel2_dir, site_code, str(year_int)) exists
                    if os.path.exists(os.path.join(sentinel2_dir, site_code, str(year_int))):
                        for file in os.listdir(os.path.join(sentinel2_dir, site_code, str(year_int))):
                            if os.path.isdir(os.path.join(sentinel2_dir, site_code, str(year_int), file)):
                                cloud_ratio = cal_cloud_covering_ratio(
                                    os.path.join(sentinel2_dir, site_code, str(year_int), file))
                                # print('Cloud ratio for %s: ' %file, cloud_ratio)
                                if cloud_ratio < cloud_ratio_threshold:
                                    sentinel2_list.append(file)
                    year_increment += 1

            if row_tocr['TOC_R'] > 0:

                year_str = row_tocr['Datetime'].split('-')[0]

                # file in sentinel2_site_dir has YYYYMMDD after the second underscore, find the file with the date closest to the datetime in tocr_data

                # Convert the date string to a datetime object
                row_date = datetime.strptime(row_tocr['Datetime'].replace('-', ''), '%Y%m%d')
                closest_file = find_closest_date_file(row_date, sentinel2_list)

                sentinel2_site_dir = os.path.join(sentinel2_dir, site_code, closest_file.split('_')[2][0:4])

                if closest_file:
                    print('Upscale datetime using Sentinel2 data: ', row_tocr['Datetime'],
                        os.path.join(sentinel2_site_dir, closest_file))
                else:
                    print('No matching file found for', row_tocr['Datetime'])

                (CGLS_grid, corrected_tocr_CGLS_resolution) = tocr_correction(os.path.join(sentinel2_site_dir, closest_file), height_tower, height_canopy, row_tocr['TOC_R'], lat, lon, OUTPUT_site_dir, row_tocr['Datetime'])

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