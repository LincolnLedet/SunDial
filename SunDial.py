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
import laspy
from pytz import timezone
from pyproj import CRS, Transformer
from scipy.ndimage import generic_filter
from pyntcloud import PyntCloud
import pandas as pd
import open3d as o3d



def view_ply(ply_path):
    # Load the point cloud
    pcd = o3d.io.read_point_cloud(ply_path)
    
    # Print info
    print(pcd)
    
    # Launch interactive viewer (mouse: orbit/zoom/pan)
    o3d.visualization.draw_geometries([pcd])

def get_now() -> datetime:
    return datetime.now(timezone('US/Eastern'))

def get_sun_alt(latitude: float, longitude: float, now: datetime) -> float:
    return pysolar.solar.get_altitude(latitude_deg=latitude, longitude_deg=longitude, when=now)

def get_sun_azu(latitude: float, longitude: float, now: datetime) -> float:
    return pysolar.solar.get_azimuth(latitude_deg=latitude, longitude_deg=longitude, when=now)

def laz_to_ply(laz_path, ply_path):
    las = laspy.read(laz_path)

    # Debug check
    print(f"x type: {type(las.x)}, shape: {np.shape(las.x)}")
    print(f"y type: {type(las.y)}, shape: {np.shape(las.y)}")
    print(f"z type: {type(las.z)}, shape: {np.shape(las.z)}")

    x = np.asarray(las.x)
    y = np.asarray(las.y)
    z = np.asarray(las.z)

    if x.ndim == 0 or y.ndim == 0 or z.ndim == 0:
        raise ValueError("Coordinate arrays are scalars. The LAZ file may not be loading properly. Check if `laspy[lazrs]` is installed.")

    data = pd.DataFrame({'x': x, 'y': y, 'z': z})
    cloud = PyntCloud(data)

    os.makedirs(os.path.dirname(ply_path), exist_ok=True)

    cloud.to_file(ply_path)
    print(f"PLY file saved to {ply_path}")

def extract_buildings_to_ply(laz_path, out_ply_path):
    import laspy
    import pandas as pd
    import numpy as np
    from pyntcloud import PyntCloud
    import os

    las = laspy.read(laz_path)

    building_mask = las.classification == 6

    if np.count_nonzero(building_mask) == 0:
        raise ValueError(f"No building points found in {laz_path}")

    x = np.asarray(las.x[building_mask])
    y = np.asarray(las.y[building_mask])
    z = np.asarray(las.z[building_mask])

    data = pd.DataFrame({'x': x, 'y': y, 'z': z})
    cloud = PyntCloud(data)

    os.makedirs(os.path.dirname(out_ply_path), exist_ok=True)
    cloud.to_file(out_ply_path)
    print(f"[+] Saved building points to {out_ply_path}")

def alpha_shape_buildings(ply_path, alpha=6.0, voxel_size=0.5):
    # Load the filtered building point cloud
    pcd = o3d.io.read_point_cloud(ply_path)
    print(f"[+] Loaded {len(pcd.points)} building points")

    # Downsample slightly to reduce noise & avoid Qhull errors
    pcd = pcd.voxel_down_sample(voxel_size=voxel_size)
    print(f"[+] Downsampled to {len(pcd.points)} points")

    if len(pcd.points) < 4:
        raise ValueError("Not enough points for surface reconstruction.")

    # Estimate normals (needed for visualization)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=1.0, max_nn=30))

    print(f"[*] Running alpha shape surface generation with alpha={alpha}")
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
    mesh.compute_vertex_normals()

    # Show result
    o3d.visualization.draw_geometries([mesh], window_name="Buildings Mesh")

    # Optionally save the mesh
    o3d.io.write_triangle_mesh("3JsFiles/building_mesh.obj", mesh)
    print("[+] Mesh saved to 3JsFiles/building_mesh.obj")

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
def fill_nodata(dem, nodata_val):
    mask = (dem == nodata_val) | np.isnan(dem)

    def mean_filter(values):
        center = values[len(values) // 2]
        if center == nodata_val or np.isnan(center):
            vals = values[(values != nodata_val) & ~np.isnan(values)]
            return np.mean(vals) if len(vals) > 0 else nodata_val
        return center

    filled = generic_filter(dem, mean_filter, size=13)
    return filled


def PlotDem(tifPath, ground = False):
    with rasterio.open(tifPath, 'r+') as src:
        dem = src.read(1)

        if ground is True:
            nodata = src.nodata if src.nodata is not None else -9999
            print("Filling missing DEM values...")
            filled_dem = fill_nodata(dem, nodata)
            src.write(filled_dem, 1)

    # Use rioxarray and earthpy to plot
    print("Plotting DEM...")
    idw_dem = rxr.open_rasterio(tifPath, masked=True)
    ep.plot_bands(
        arr=idw_dem,
        cmap="RdYlGn",
        figsize=(8, 8),
        title="LiDAR DEM (Filled)"
    )
    plt.show()


#Creates DEM .tif files
def RenderLidar(laz_file_path, out_dir, out_name,excludes, returns="all", ):
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
                if (val > 245):
                    #print(val)
                    dem[row, col] = 500

        # Write our modified array back into band #1
        src.write(dem, 1)


def main():
    latitude = 33.94978248553038
    longitude = -83.37334782237981
    cwd =  os.getcwd()
    interpolated_directory = os.path.join(cwd, "output")
    mapPath = os.path.join(cwd, fetch_lidar_file(latitude,longitude))

    extract_buildings_to_ply(mapPath, os.path.join(cwd, "3JsFiles", "pointcloud.ply"))
    alpha_shape_buildings("3JsFiles/pointcloud.ply")

    print(mapPath)
    groundMap = os.path.join(cwd, "output\groundMap.tif")
    BuildingMap = os.path.join(cwd, "output\BuildingMap.tif")

    #RenderLidar(laz_file, interpolated_directory, "TallCaster.tif", "0,1,2,3,4,7,8,9,12,13,14", returns="all")
    #RenderLidar(laz_file, interpolated_directory, "ShortCaster.tif", "0,1,2,6,7,8,9,10,11,12,13,14", returns="all")
    #print_laspy_info(laz_file)
    #RenderLidar(mapPath, interpolated_directory, "groundMap.tif","3,4,5,6", returns="all")
    #PlotDem(groundMap, ground = True)
    #RenderLidar(mapPath, interpolated_directory, "BuildingMap.tif","1,2,3,4,5", returns="all")
    #PlotDem(BuildingMap, ground = False)
   # RenderLidar(mapPath, interpolated_directory, "map.tif","0", returns="all")
    #PlotDem(map)
    #rasterioTest(ShortCaster)

    #rasterioTest(SurfacePath)
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
    