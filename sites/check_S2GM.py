#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    create_AOIs.py
# @Author:      Dr. Jean-Sébastien Carrière
# @Email:       jean-sebastien.carriere@acri-st.fr
# @Time:        26/09/2024 17:56

# This script is used to check the input S2GM products.

import os
import csv
import glob
import time
import shutil
import argparse

from IPython import embed

S2GM_path = "/mount/internal/work-st/projects/jrc-066/1953-s2gm/SIAC-Albedo/Data/S2GM/GBOV"
S2GM_QLs = "/mount/internal/work-st/projects/jrc-066/1953-s2gm/SIAC-Albedo/Data/S2GM/GBOV/QLs"
if not os.path.isdir(S2GM_QLs): os.mkdir(S2GM_QLs)
else: 
    shutil.rmtree(S2GM_QLs)
    time.sleep(2)
    os.mkdir(S2GM_QLs)

parser = argparse.ArgumentParser(description="Select the Site you want to check")
parser.add_argument('-s', '--site', type=str, default='', help="Specific site")

args = parser.parse_args()
Site = args.site
Year = "2022"

S2GM_products = glob.glob(os.path.join(S2GM_path, Site, Year, "*"))
Empty_products = []
for product in S2GM_products:
    QLs = glob.glob(os.path.join(product, 'QuickLookRGB*'))
    if len(QLs) == 6:
        for QL in QLs:
            shutil.copy(QL, os.path.join(S2GM_QLs, QL.split(f"_{Site}_")[-1].split('_')[0] + '_' + os.path.basename(QL)))
    else:
        Empty_products.append(product.split('_')[-3])
        shutil.rmtree(product)


print(Empty_products)
embed()