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


def fetch_lidar_file(latidtue, longitude): #pulls .laz file from usgs
    url = "https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/Projects/GA_Statewide_2018_B18_DRRA/GA_Statewide_B4_2018/LAZ/USGS_LPC_GA_Statewide_2018_B18_DRRA_e1157n1284.laz"
    file_name = url.split("/")[-1]
    output_path = os.path.join("LAZ", file_name)
    # Create the output directory if it doesn't exist
    os.makedirs("LAZ", exist_ok=True)

    print(f"Starting download: {file_name}")
    with requests.get(url, stream=True, verify=False) as response:
        response.raise_for_status()
        with open(output_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
    print(f"Download complete: {output_path}")

    return

# prints raw data from LAS file. Note that a LAZ file is a compressed LAS file
def print_laspy_info(laz_file_path):
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
        resolution="1",
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
                #if (val < 0):
                    #print(val)
                #dem[row, col] = 0

        # Write our modified array back into band #1
        src.write(dem, 1)
    
    


def main():
    #fetch_lidar_file(1, 1)
    interpolated_directory = r"C:\Users\lll81910\Desktop\Coding Projects\SunDial\output"
    laz_file = r"C:\Users\lll81910\Desktop\Coding Projects\SunDial\LAZ\USGS_LPC_GA_Statewide_2018_B18_DRRA_e1157n1283.laz"
    TallCaster = r"C:\Users\lll81910\Desktop\Coding Projects\SunDial\output\TallCaster.tif"
    ShortCaster = r"C:\Users\lll81910\Desktop\Coding Projects\SunDial\output\ShortCaster.tif"

    SurfacePath = r"C:\Users\lll81910\Desktop\Coding Projects\SunDial\output\Surface.tif"



    #RenderLidar(laz_file, interpolated_directory, "TallCaster.tif", "0,1,2,3,4,7,8,9,12,13,14", returns="all")
    RenderLidar(laz_file, interpolated_directory, "ShortCaster.tif", "0,1,2,6,7,8,9,10,11,12,13,14", returns="all")
    #RenderLidar(laz_file, interpolated_directory, "Surface.tif","0,1,3,4,5,6,7,8,9,12,13,14", returns="all")
    rasterioTest(TallCaster)
    rasterioTest(ShortCaster)
    rasterioTest(SurfacePath)

   # rasterioTest(SurfacePath)
    tifEditTest(TallCaster)



    output_directory = r"C:\Users\lll81910\Desktop\Coding Projects\SunDial\output"
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

