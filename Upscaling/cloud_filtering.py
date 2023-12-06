#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    cloud_filtering.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        06/12/2023 15:47

from osgeo import gdal
import numpy as np
import os

def cal_cloud_covering_ratio(sentinel2_dir):

    sentinel2_dhr_dir = os.path.join(sentinel2_dir, 'tile_0', 'albedo')
    for dhr_file in os.listdir(sentinel2_dhr_dir):
        if dhr_file.endswith('B02_UCL_dhr.jp2'):
            dhr_b02 = gdal.Open(os.path.join(sentinel2_dhr_dir, dhr_file))

    # get the center index of the image.
    center_index_x = int(dhr_b02.RasterXSize / 2)
    center_index_y = int(dhr_b02.RasterYSize / 2)

    # get the subarray from plus minus 200 pixels from the center.
    subarray = dhr_b02.ReadAsArray(center_index_x - 200, center_index_y - 200, 400, 400)
    cloud_ratio = np.size(subarray[subarray < 0]) / np.size(subarray)
