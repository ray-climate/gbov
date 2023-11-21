#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    create_AOIs.py
# @Author:      Dr. Rui Song
# @Email:       rui.song@physics.ox.ac.uk
# @Time:        21/11/2023 14:45

# This script is used to create AOIs for a given latitude and longitude,
# generate a geojson file 20km * 20km centered at the given lat and lon

import csv
import json
import os
from geopy.distance import great_circle

def create_geojson(lat, lon, size_km, filename, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Function to calculate the offset in degrees for a given distance in km
    def calculate_offset(lat, lon, distance_km, direction):
        directions = {
            'north': 0,
            'east': 90,
            'south': 180,
            'west': 270
        }
        bearing = directions[direction]
        destination_point = great_circle(kilometers=distance_km).destination((lat, lon), bearing)
        if direction in ['north', 'south']:
            return destination_point.latitude - lat
        else:
            return destination_point.longitude - lon

    half_size = size_km / 2

    # Calculate bounds
    north_offset = calculate_offset(lat, lon, half_size, 'north')
    south_offset = calculate_offset(lat, lon, half_size, 'south')
    east_offset = calculate_offset(lat, lon, half_size, 'east')
    west_offset = calculate_offset(lat, lon, half_size, 'west')

    # Define corners
    northwest = [lon + west_offset, lat + north_offset]
    northeast = [lon + east_offset, lat + north_offset]
    southeast = [lon + east_offset, lat + south_offset]
    southwest = [lon + west_offset, lat + south_offset]

    # Create GeoJSON Feature
    feature = {
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "Polygon",
            "coordinates": [[northwest, northeast, southeast, southwest, northwest]]
        }
    }

    # Create GeoJSON object
    geojson_obj = {
        "type": "FeatureCollection",
        "features": [feature]
    }

    # Write to file in the specified output directory
    output_path = os.path.join(output_dir, f'{filename}.geojson')
    with open(output_path, 'w') as file:
        json.dump(geojson_obj, file, indent=4)

    print(f"GeoJSON file '{output_path}' created successfully.")

# Read the CSV and process each row
def process_csv(file_path, size_km, output_dir):
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            location, lat, lon = row
            lat, lon = float(lat), float(lon)  # Convert latitude and longitude to float
            create_geojson(lat, lon, size_km, location, output_dir)

# Example usage
output_directory = 'geojson_output'  # Define your output directory
process_csv('../sites/geolocation.csv', 20, output_directory)  # Adjust the size_km as needed
