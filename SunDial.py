#Code by Lincoln Ledet

import whitebox
import whitebox_workflows as wbw
import geemap
import rioxarray as rxr
import rasterio
import earthpy as et
import earthpy.plot as ep
import matplotlib.pyplot as plt
import geopandas as gpd
from rasterio.plot import plotting_extent
from pyproj import Transformer

#All local uga lidar scans can be found here 
#https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/Projects/GA_Statewide_2018_B18_DRRA/GA_Statewide_B4_2018/LAZ/
# These file names use EPSG:6350 / Conus Albers naming format
# Things like google maps use 4326 WGS 84 naming format
# here is a calculator https://epsg.io/transform#s_srs=4326&t_srs=6350&x=-83.3487450&y=33.9212080


# Why does Whitebox tools fucking rock you ask
# its a free open-source, cross-platfrom geospatial anlysis library written in rust.
# THIS MF DOES RASTERS, VECTORS, WATERSHEDS. and most importantly LIDAR POINT CLOUD OPS!
# C:\Users\linco\Desktop\Ledet_Code\SunDial\output


wbt = whitebox.WhiteboxTools() 
wbt.set_working_dir(r"C:\Users\linco\Desktop\Ledet_Code\SunDial") #set working direectory

wbt.lidar_idw_interpolation( # this takes raw lidar points and interpolates a continuse raster using inverse-distance weighting 
    i="USGS_LPC_GA_Statewide_2018_B18_DRRA_e1157n1283.laz", #this is where the LAZ file lives
    output=r"C:\Users\linco\Desktop\Ledet_Code\SunDial\output\interpolated.tif",# write path
    parameter="elevation", #intensity, number_of_returns 
    returns="all", # ground(last) vs. vegetation(first) differentiation
    resolution=".5",# grid cell size in metres (1 m x 1 m pixels).
    exclude_cls="1,2,3,4,5,7,8,9,10,12,13" # what classes we want to exclude (see list at bottom of file). Note this if all forms of data are excluded the tif file wont render
)

geemap.add_crs(r"C:\Users\linco\Desktop\Ledet_Code\SunDial\output\interpolated.tif",
               epsg=32613) # This CRS is based on one of the shapefiles within the downloaded data

idw_dem = rxr.open_rasterio(r"C:\Users\linco\Desktop\Ledet_Code\SunDial\output\interpolated.tif",
                            masked=True)

ep.plot_bands(
    arr=idw_dem,
    cmap="RdYlGn",
    figsize=(8, 8),
    title="Lidar DEM using IDW"
)

plt.show()
print("🚀 All geospatial tools imported successfully!")



# ASPRS LAS classification table (most common codes you’ll see):
# Code  Class                       Notes
# 0     Created, never classified   Default “no idea”
# 1     Unclassified                Points yet to be classified
# 2     Ground                      Bare-earth returns
# 3     Low vegetation              Shrubs, grass (<0.5 m tall)
# 4     Medium vegetation           Brush (0.5 m–2 m tall)
# 5     High vegetation             Trees (>2 m tall)
# 6     Building                    Roof tops, structures
# 7     Low point (key point)       Tie-points for breaklines
# 8     Model key-point             Critical edge points
# 9     Water                       Lakes, rivers
# 10    Rail                        Railroad surfaces
# 11    Road surface                Paved roads
# 12    Bridge deck                 Bridge tops
# 13    High noise                  Bad returns, birds, etc.
