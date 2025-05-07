Overview
The goal of this project is to create a publicly accessible, open-source API or Python library that can take in the input longitude, latitude, and time of day, and predict if a point is in shade based on publicly available LiDAR scans. This project leverages powerful geospatial libraries for LiDAR data processing and analysis.


Requirements
To run this project on Windows, the following packages are required:
    whitebox - For LiDAR point cloud processing and terrain analysis.
    whitebox_workflows - High-level geospatial analysis tools.
    geemap - Simplified Google Earth Engine Python interface.
    rioxarray - Raster file processing with support for CRS transformations.
    rasterio - Reading and writing geospatial raster data.
    earthpy - Tools for spatial data analysis and visualization.
    matplotlib - Plotting library for data visualization.
    geopandas - Geospatial data processing with pandas-like syntax.
    pyproj - Coordinate transformation and CRS conversions

Installation
First, create a new Conda environment to isolate the project dependencies:
conda create -n shade-prediction python=3.11
conda activate shade-prediction
conda install -c conda-forge whitebox whitebox_workflows geemap rioxarray rasterio earthpy matplotlib geopandas pyproj

Lidar scans of UGA can be found here 
https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/Projects/GA_Statewide_2018_B18_DRRA/GA_Statewide_B4_2018/LAZ/

All local UGA LiDAR scans can be found here. These files use EPSG:6350 (CONUS Albers) as the coordinate system.  
These file names use EPSG:6350 / Conus Albers naming format
Things like google maps use 4326 WGS 84 naming format
here is a calculator https://epsg.io/transform#s_srs=4326&t_srs=6350&x=-83.3487450&y=33.9212080


