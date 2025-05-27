#Code by Lincoln Ledet

import os
import requests
import whitebox
import geemap
import rioxarray as rxr
import earthpy.plot as ep
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import pysolar
import datetime
from pytz import timezone
from pyproj import CRS, Transformer


def get_now() -> datetime:
    return datetime.now(timezone('US/Eastern'))

def get_sun_alt(latitude: float, longitude: float, now: datetime) -> float:
    return pysolar.solar.get_altitude(latitude_deg=latitude, longitude_deg=longitude, when=now)

def get_sun_azu(latitude: float, longitude: float, now: datetime) -> float:
    return pysolar.solar.get_azimuth(latitude_deg=latitude, longitude_deg=longitude, when=now)

def fetch_lidar_file(latitude, longitude): #return file name
    #dealing with location format (coverting latitude and longitude to Conus Albers)
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

# Plots the DEM
def rasterioTest(tifPath):
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


#Creates DEM .tif files
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
        resolution="2",
        exclude_cls= excludes
    )


#edits raw DEM grid
def tifEditTest(TifPath):
    with rasterio.open(TifPath, mode='r+') as src:
    # Read the first (and usually only) band into a NumPy array
        dem = src.read(1)

        # Loop over every row and column
        print (src.height,src.width)
        for row in range(src.height):
            for col in range(src.width):
                val = dem[row, col]
                if (val > 245):
                    #print(val)
                    dem[row, col] = 500

        # Write our modified array back into band #1
        src.write(dem, 1)


def main():
    latitude = 33.9224066505482
    longitude = -83.35140814950346

    cwd =  os.getcwd()
    interpolated_directory = os.path.join(cwd, "output")
    mapPath = os.path.join(cwd, fetch_lidar_file(latitude,longitude))
    print(mapPath)

    map = os.path.join(cwd, "output\map.tif")

    #RenderLidar(laz_file, interpolated_directory, "TallCaster.tif", "0,1,2,3,4,7,8,9,12,13,14", returns="all")
    #RenderLidar(laz_file, interpolated_directory, "ShortCaster.tif", "0,1,2,6,7,8,9,10,11,12,13,14", returns="all")
    #print_laspy_info(laz_file)
    RenderLidar(mapPath, interpolated_directory, "map.tif","0", returns="all")
    #rasterioTest(TallCaster)
    #rasterioTest(ShortCaster)
    rasterioTest(map)

   # rasterioTest(SurfacePath)
    #plot_lidar_file(laz_file, output_directory, crs=6350, resolution=0.5)


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
# 12: Reserved
# 13: Wire - Guard (Shield)
# 14: Wire - Conductor (Phase)
# 15: Transmission Tower
# 16: Wire-Structure Connector (Insulator)
# 17: Bridge Deck
# 18: High Noise  


# TO DO
# Add comments and use os.path to make code more accessible. 
# Finish Method to download usgs laz file based on lat and long
# Preprocess laz and las 
    