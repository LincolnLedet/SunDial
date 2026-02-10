#Code by Lincoln Ledet
import os
import requests
import whitebox
import geemap
import rioxarray as rxr
import earthpy.plot as ep
import matplotlib.pyplot as plt
import ssl
import laspy
import numpy as np
import rasterio
from pyproj import CRS, Transformer
from scipy.interpolate import griddata
from datetime import datetime
import pandas as pd
from pvlib import solarposition
import math
import pytz
from numba import njit


def fetch_lidar_file(latitude, longitude): 
    """
    Downloads Compressed Lidar Data (laz) from usgs

    Parameters:
    latitude (float)
    latitude (float)

    Returns:
    String: Download location

    Raises:
    """
    alb   = CRS.from_epsg(6350)
    wgs84 = CRS.from_epsg(4326)
    to_alb = Transformer.from_crs(wgs84, alb, always_xy=True)
    albx, alby = to_alb.transform(longitude, latitude)
    albx = int(albx // 1000) 
    alby = int(alby // 1000) # 
    #url construction
    url = "https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/Projects/GA_Statewide_2018_B18_DRRA/GA_Statewide_B4_2018/LAZ/USGS_LPC_GA_Statewide_2018_B18_DRRA_e" + str(albx) + "n" + str(alby) + ".laz"
    print("fetching laz file from " + url)
    file_name = url.split("/")[-1]
    output_path = os.path.join("LAZ", file_name)
    os.makedirs("LAZ", exist_ok=True)    # Create the output directory if it doesn't exist

    print(f"Starting download: {file_name}")
    with requests.get(url, stream=True, verify=False) as response:
        response.raise_for_status()
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
    print(f"Download complete: {output_path}")
    return (output_path)



def print_laspy_info(laz_file_path):
    """
    Prints laz file data

    Parameters:
    laz_file_path (string)
    latitude (float)

    """
    with laspy.open(laz_file_path) as fh:
        print('Points from Header:', fh.header.point_count)
        las = fh.read()

        print(las)
        print('Points from data:', len(las.points))
        ground_pts = las.classification == 2
        bins, counts = np.unique(las.return_number[ground_pts], return_counts=True)
        print('Ground Point Return Number distribution:')
        for r,c in zip(bins,counts):
            print('    {}:{}'.format(r,c))

    return

# Plots the DEM
def rasterioTest(tifPath):
    """
    Creates a DEM plot from a tif file

    Parameters:
    tifPath (String)

    """
    src= rasterio.open(tifPath)
    dem_data = src.read(1)
    print(src.meta)

    # Plot the resulting raster
    print("Plotting DEM...")
    idw_dem = rxr.open_rasterio(tifPath, masked=True)
    ep.plot_bands(
        arr=idw_dem,
        cmap="RdYlGn",
        figsize=(8, 8),
        title="LiDAR DEM using IDW"
    )
    plt.show()
    return


def RenderLidar(laz_file_path, out_dir, out_name,excludes, returns="all"):
    """
    runs lidar_idw_interpolation:
      - returns: "all" or "first"
      - keep_classes: list of ASPRS codes to INCLUDE (others get excluded)
    """
    wbt = whitebox.WhiteboxTools()
    wbt.set_working_dir(out_dir)
# ONLY exclude if keep_classes is explicitly given
    
    if not os.path.exists(laz_file_path):
        raise FileNotFoundError(f"LAZ file not found: {laz_file_path}")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    wbt.lidar_idw_interpolation(
        i=laz_file_path,
        output=os.path.join(out_dir, out_name),
        parameter="elevation",
        returns=returns,
        resolution=".25",
        exclude_cls= excludes
    )


def fill_dem_holes(input_tif: str, output_tif: str):
    """
    runs fill_missing_data:
      - Creates a new tif.with filled in holes
    """
    wbt = whitebox.WhiteboxTools()
    wbt.set_working_dir(os.path.dirname(input_tif))
    wbt.fill_missing_data(
        i=input_tif,
        output=output_tif,
        filter=61,      # kernel size, must be odd
        weight=2.0      # how strongly closer points are weighted
    )


def fill_dem_gaps(input_tif, output_tif):
    """
    runs fill_missing_data:
      - Creates a new tif.with filled in holes
    """
    with rasterio.open(input_tif) as src:
        dem = src.read(1)
        profile = src.profile.copy()
        nodata_val = src.nodata if src.nodata is not None else -9999

    # Replace NoData values with np.nan
    dem = np.where(dem == nodata_val, np.nan, dem)

    # Generate coordinate grid
    rows, cols = np.indices(dem.shape)
    coords = np.column_stack((rows.ravel(), cols.ravel()))

    # Get valid and missing point positions
    valid_mask = ~np.isnan(dem)
    missing_mask = np.isnan(dem)

    valid_coords = coords[valid_mask.ravel()]
    valid_values = dem[valid_mask]

    missing_coords = coords[missing_mask.ravel()]

    # Run interpolation (nearest fills large gaps best)
    filled_values = griddata(
        valid_coords,
        valid_values,
        missing_coords,
        method='nearest'  # 'linear' or 'cubic' are smoother but need dense data
    )

    # Inject filled values
    dem[missing_mask] = filled_values

    # Write result
    with rasterio.open(output_tif, 'w', **profile) as dst:
        dst.write(dem, 1)



@njit
def _cast_shadows(dem, height, width, dx, dy, tan_el, resolution, max_distance):
    steps = int(max_distance / resolution)  # how many pixels to check
    for row in range(height - 200):
        for col in range(width - 200):
            z0 = dem[row, col]  # starting elevation
            for k in range(1, steps):
                # project ray from current point
                xi = int(round(col + k * dx))
                yi = int(round(row + k * dy))

                #check bounds
                if 0 <= yi < height and 0 <= xi < width:
                    dist_m = k * resolution  # actual distance in meters
                    z_ray = z0 + dist_m * tan_el  # ray height at that distance
                    if dem[yi, xi] > z_ray:  # terrain blocks the sun
                        dem[row, col] = 160  # mark shadow
                        break
                else:
                    break  # ray has left the grid
    return dem

def ground_level_shadow_cast(TifPath, latitude, longitude):

    # Get current time in Eastern
    eastern = pytz.timezone('America/New_York')
    now = pd.Timestamp(datetime.now(eastern))
    solpos = solarposition.get_solarposition(time=now, latitude=latitude, longitude=longitude)

    # Extract elevation and azimuth
    elevation = solpos['elevation'].iloc[0]
    azimuth = solpos['azimuth'].iloc[0]

    print(f"Sun Elevation: {elevation:.2f}°")
    print(f"Sun Azimuth: {azimuth:.2f}°")

    azrad = np.deg2rad(azimuth)
    elrad = np.deg2rad(elevation)
    dx = np.sin(azrad)  # how much to move east per step and not 0 degree is north
    dy = -np.cos(azrad)  # negated because raster rows increase southward

    tan_el = np.tan(elrad) # height gain per unit horizontal distance

    with rasterio.open(TifPath, mode='r+') as src:
        resolution = src.res[0]  # pixel size in meters from the GeoTIFF
        dem = src.read(1)
        dem = _cast_shadows(dem, src.height, src.width, dx, dy, tan_el, resolution, 25.0)
        src.write(dem, 1)
    return

    
    
def main():
    latitude = 33.96280751837907
    longitude = -83.38712753678465

    cwd =  os.getcwd()
    interpolated_directory = os.path.join(cwd, "output")
    laz_file = os.path.join(cwd, fetch_lidar_file(latitude, longitude))

    #Render Dem model from lidar point cloud
    RenderLidar(laz_file, interpolated_directory, "Ground.tif", "0,1,3,4,5", returns="all")
    
    #init file paths
    Ground_raw = os.path.join(interpolated_directory, "Ground.tif")
    Ground_filled = os.path.join(interpolated_directory, "Ground_filled.tif")
    Ground_filled_again = os.path.join(interpolated_directory, "Ground_filled_again.tif")

    #Fills small holes in lidar data (accurate)
    fill_dem_holes(Ground_raw, Ground_filled)

    #Fills big gaps in lidar data (less accurate)
    fill_dem_gaps(Ground_filled, Ground_filled_again)

    #displays non 
    rasterioTest(Ground_filled_again)

    #lets try to cast some shadows!    
    ground_level_shadow_cast(Ground_filled_again, latitude, longitude)

    #display dem file with shadows casted
    rasterioTest(Ground_filled_again)



if __name__ == "__main__":
    main()

# ASPRS LAS classification mapping:
# 0:  Never classified
# 1:  Unassigned
# 2:  Ground
# 3:  Low vegetation
# 4:  Medium vegetation
# 5:  High vegetation
# 6:  Building
# 7:  Low point (noise)
# 8:  Reserved
# 9:  Water
# 10: Rail
# 11: Road surface


    