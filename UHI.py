#!/usr/bin/env python
# coding: utf-8

# Import necessary libraries
import ee
import geemap
import geopandas as gpd

# Initialize the Earth Engine module
ee.Authenticate()
ee.Initialize()

#========================

# AOI

# Load the AOI GeoJSON file as a GeoDataFrame
geojson_url = 'https://drive.google.com/uc?id=1O9nrrGjMIX3naiHPnsmXS2X-QWQfWdIu'
gdf = gpd.read_file(geojson_url)

# Convert the GeoDataFrame to an Earth Engine FeatureCollection
aoi = geemap.geopandas_to_ee(gdf)

#========================

# Cloud mask function for Sentinel-2 Harmonized
def mask_cloud_s2_harmonized(image):
    # Sentinel-2 QA60 band: clouds and cirrus
    qa = image.select('QA60')
    cloud_mask = qa.bitwiseAnd(1 << 10).eq(0).And(qa.bitwiseAnd(1 << 11).eq(0))  # Mask clouds
    return image.updateMask(cloud_mask).divide(10000).select(
        ['B4', 'B8', 'B11'],  # Red, NIR, SWIR bands
        ['B4', 'B8', 'B11']   # Rename bands for consistency
    )

# Sentinel-2 Harmonized Image Collection
sentinel2_harmonized_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')

# Define specific imagery dates
dates = ['2024-08-30', '2024-08-25', '2024-08-20', '2024-08-15', '2024-07-11', 
         '2024-07-06', '2024-06-26', '2024-06-16', '2024-06-11']

# Convert dates to a list of Earth Engine-compatible milliseconds
date_millis = [ee.Date(date).millis() for date in dates]

# Filter images based on date list
date_filter = ee.Filter.inList('system:time_start', date_millis)

filtered_collection = (
    sentinel2_harmonized_collection
    .filter(date_filter)
    .filterBounds(aoi)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(mask_cloud_s2_harmonized)
)

#========================

# Function to calculate NDVI
def calculate_ndvi(image):
    ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')  # NIR and Red bands
    return image.addBands(ndvi)

# Apply NDVI calculation to each image in the collection
ndvi_collection = filtered_collection.map(calculate_ndvi)

#========================

# Brightness Temperature from Sentinel-2 Harmonized
def calculate_lst(image):
    # Use NDVI from processed bands
    ndvi = image.select('NDVI')
    
    # Calculate PV (Proportion of Vegetation)
    pv = ndvi.subtract(0.2).divide(0.5).clamp(0, 1).pow(2).rename('PV')  # PV between 0 and 1
    
    # Calculate Emissivity
    em = pv.multiply(0.004).add(0.986).rename('EM')  # Simplified formula for emissivity
    
    # Convert SWIR band (B11) to Brightness Temperature (assuming TOA radiance approximation)
    brightness_temp = image.select('B11').multiply(1e4).rename('BT')  # Scale to TOA radiance
    
    # Calculate LST in Celsius
    lst = brightness_temp.expression(
        '(BT / (1 + (0.00115 * (BT / 1.438)) * log(em))) - 273.15',
        {'BT': brightness_temp, 'em': em}
    ).rename('LST')
    
    return lst

# Apply LST calculation to each image in the collection
lst_collection = ndvi_collection.map(calculate_lst)

#========================

# Calculate the median LST over the entire period
median_lst = lst_collection.median().clip(aoi)

# Initialize a map centered around the AOI
Map = geemap.Map()
Map.centerObject(aoi, zoom=10)  # Adjust zoom level as needed

# Visualization parameters for LST
lst_vis_params = {
    'min': 20,
    'max': 48,
    'palette': ['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
}

# Add the AOI, mean LST, and other layers to the map
Map.addLayer(aoi, {}, "AOI")
Map.addLayer(median_lst, lst_vis_params, 'Median LST AOI')
Map

#========================

# Normalized UHI
lst_mean = median_lst.reduceRegion(
    reducer=ee.Reducer.mean(),
    geometry=aoi,
    scale=10,
    maxPixels=1e9
).values().get(0)

lst_std = median_lst.reduceRegion(
    reducer=ee.Reducer.stdDev(),
    geometry=aoi,
    scale=10,
    maxPixels=1e9
).values().get(0)

# Convert lst_mean and lst_std to ee.Image objects
lst_mean_img = ee.Image.constant(lst_mean)
lst_std_img = ee.Image.constant(lst_std)

print('Mean LST in AOI:', lst_mean.getInfo())
print('STD LST in AOI:', lst_std.getInfo())

# Calculate UHI by subtracting and dividing as needed
uhi = median_lst.subtract(lst_mean_img).divide(lst_std_img).rename('UHI')

uhi_clipped = uhi.clip(aoi)

# Define visualization parameters
uhi_clipped_vis = {
    'min': -4,
    'max': 4,
    'palette': ['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
}

# Add UHI layer to the map
Map.addLayer(uhi_clipped, uhi_clipped_vis, 'UHI AOI')

# Urban Thermal Field Variance Index (UTFVI)
utfvi = median_lst.subtract(lst_mean_img).divide(median_lst).rename('UTFVI')

utfvi_clipped = utfvi.clip(aoi)

utfvi_clipped_vis = {
    'min': -1,
    'max': 1,
    'palette': ['313695', '74add1', 'fed976', 'feb24c', 'fd8d3c', 'fc4e2a', 'e31a1c', 'b10026']
}
Map.addLayer(utfvi_clipped, utfvi_clipped_vis, 'UTFVI AOI')
Map


# In[ ]:




