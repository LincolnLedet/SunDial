#Code by Lincoln Ledet

import os
import whitebox
import geemap
import rioxarray as rxr
import earthpy.plot as ep
import matplotlib.pyplot as plt


def process_lidar_file(laz_file_path, output_dir, crs=6350, resolution=0.5, exclude_classes="1,2,3,4,5,7,8,9,10,12,13"):
    # Validate file paths
    if not os.path.exists(laz_file_path):
        raise FileNotFoundError(f"LAZ file not found: {laz_file_path}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Set up Whitebox tools
    wbt = whitebox.WhiteboxTools()
    wbt.set_working_dir(output_dir)

    # Generate output file path
    output_tif = os.path.join(output_dir, "interpolated.tif")

    # Run IDW interpolation
    print("Starting LiDAR interpolation...")
    wbt.lidar_idw_interpolation(
        i=laz_file_path,
        output=output_tif,
        parameter="elevation",
        returns="all",
        resolution=str(resolution),
        exclude_cls=exclude_classes
    )

    # Add CRS to the output file
    geemap.add_crs(output_tif, epsg=crs)

    # Plot the resulting raster
    print("Plotting DEM...")
    idw_dem = rxr.open_rasterio(output_tif, masked=True)
    ep.plot_bands(
        arr=idw_dem,
        cmap="RdYlGn",
        figsize=(8, 8),
        title="LiDAR DEM using IDW"
    )
    plt.show()

    print("🚀 LiDAR processing complete!")


def main():
    laz_file = r"C:\Users\linco\Desktop\Ledet_Code\SunDial\USGS_LPC_GA_Statewide_2018_B18_DRRA_e1157n1283.laz"
    output_directory = r"C:\Users\linco\Desktop\Ledet_Code\SunDial\output"
    process_lidar_file(laz_file, output_directory, crs=6350, resolution=0.5)


if __name__ == "__main__":
    main()