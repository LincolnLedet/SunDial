# Shade Prediction API

## Overview
The goal of this project is to create a publicly accessible, open-source API or Python library that can take in the input longitude, latitude, and time of day, and predict if a point is in shade based on publicly available LiDAR scans. This project leverages powerful geospatial libraries for LiDAR data processing and analysis.

---

## Requirements
To run this project on Windows, the following packages are required:

- **whitebox** - For LiDAR point cloud processing and terrain analysis.  
- **whitebox_workflows** - High-level geospatial analysis tools.  
- **geemap** - Simplified Google Earth Engine Python interface.  
- **rioxarray** - Raster file processing with support for CRS transformations.  
- **rasterio** - Reading and writing geospatial raster data.  
- **earthpy** - Tools for spatial data analysis and visualization.  
- **matplotlib** - Plotting library for data visualization.  
- **geopandas** - Geospatial data processing with pandas-like syntax.  
- **pyproj** - Coordinate transformation and CRS conversions.  

---

## Installation
First, create a new Conda environment to isolate the project dependencies:

```bash
conda create -n SunDial_env python=3.11
conda activate SunDial_env
conda install -c conda-forge whitebox whitebox_workflows geemap rioxarray rasterio earthpy matplotlib geopandas pyproj
