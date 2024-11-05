#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    launcher.py
# @Author:      Dr. Jean-Sébastien Carrière
# @Email:       jean-sebastien.carriere@acri-st.fr
# @Time:        14/10/2024 17:59

from ReferenceMeasurements import preprocess_reference_measurements, process_all_sites
from Upscaling import cal_upscaling, QCheck

def main():
    print("Launching all scripts")
    preprocess_reference_measurements.main()
    process_all_sites.main()
    cal_upscaling.main()
    QCheck.main()

if __name__ == "__main__":
    main()