# gbov

# Overview
The gbov processing chain consists of a series of Python scripts designed to process and analyse shortwave broadband raidiation data from tower albedometre and satellite data. The primary focus is on preprocessing and upscaling BRF/albedo from tower to CGLS coarse resolutions.

# Features

create_AOIs.py: This script generates a 20km by 20km subset geojson file, serving as an input area for S2GM processing.

preprocess_reference_measurements.py: This script preprocesses tower data acquired from the GbOV site RM1. It filters the data and extracts key measurements like SW_IN, SW_OUT, SW_DIF.

process_all_sites.py: The core functionality of this script involves calculating point Top-of-Canopy Reflectance (ToC-R), Direct Hemispherical Reflectance (DHR), and Bi-Hemispherical Reflectance (BHR) using filtered raw measurements.

cal_upscaling.py: Assuming the point ToC-R, DHR, and BHR values from tower measurements are accurate, this script uses S2GM calculated Bidirectional Reflectance Factor (BRF) and albedo images to calculate upscaled values at CGLS grid resolution.
