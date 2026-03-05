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
from datetime import datetime, timedelta
import pandas as pd
from pvlib import solarposition
import math
import pytz
from numba import njit, prange
import imageio
import shutil


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


def clip_laz_to_radius(laz_file_path, latitude, longitude, radius_ft=100):
    """Clips a LAZ file to a circle around the given lat/long. Returns path to clipped file."""
    radius_m = radius_ft * 0.3048

    # Transform lat/long to the LAZ file's CRS
    laz = laspy.read(laz_file_path)
    laz_crs = CRS.from_wkt(laz.header.parse_crs().to_wkt()) if laz.header.parse_crs() else CRS.from_epsg(6350)
    wgs84 = CRS.from_epsg(4326)
    to_laz = Transformer.from_crs(wgs84, laz_crs, always_xy=True)
    cx, cy = to_laz.transform(longitude, latitude)

    # Filter points within radius
    dx = laz.x - cx
    dy = laz.y - cy
    dist = np.sqrt(dx**2 + dy**2)
    mask = dist <= radius_m

    clipped = laspy.LasData(laz.header)
    clipped.points = laz.points[mask]

    clipped_path = laz_file_path.replace(".laz", "_clipped.laz")
    clipped.write(clipped_path)
    print(f"Clipped {mask.sum()} of {len(laz.points)} points within {radius_ft}ft radius")
    return clipped_path

def RenderLidar(laz_file_path, out_dir, out_name,excludes, returns="all"):
    """
    runs lidar_tin_gridding:
      - returns: "all" or "first"
      - excludes: ASPRS classification codes to exclude
    """
    wbt = whitebox.WhiteboxTools()
    wbt.set_working_dir(out_dir)

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



@njit(cache=True)
def _sample_bilinear(dem, x, y, height, width):
    """Bilinear interpolation of DEM at fractional (x, y) coordinates."""
    x0 = int(x)
    y0 = int(y)
    x1 = min(x0 + 1, width - 1)
    y1 = min(y0 + 1, height - 1)
    fx = x - x0
    fy = y - y0
    return (dem[y0, x0] * (1 - fx) * (1 - fy) +
            dem[y0, x1] * fx * (1 - fy) +
            dem[y1, x0] * (1 - fx) * fy +
            dem[y1, x1] * fx * fy)

@njit(parallel=True, cache=True)
def _cast_shadows(dem, height, width, dx, dy, tan_el, resolution, max_distance):
    steps = int(max_distance / resolution)  # how many pixels to check
    shadow = np.zeros((height, width), dtype=np.uint8)
    for row in prange(height):
        for col in range(width):
            z0 = dem[row, col]
            for k in range(1, steps):
                x = col + k * dx
                y = row + k * dy

                if 1 <= y < height - 1 and 1 <= x < width - 1:
                    dist_m = k * resolution
                    z_ray = z0 + dist_m * tan_el
                    terrain_h = _sample_bilinear(dem, x, y, height, width)
                    if terrain_h > z_ray:
                        shadow[row, col] = 1
                        break
                else:
                    break
    return shadow

def ground_level_shadow_cast(TifPath, latitude, longitude, time_dt=None):

    eastern = pytz.timezone('America/New_York')
    if time_dt is None:
        now_dt = datetime.now(eastern)
        if now_dt.hour >= 18 and now_dt.minute >= 30 or now_dt.hour >= 19:
            now_dt = now_dt.replace(hour=12, minute=0, second=0)
    else:
        now_dt = time_dt
    now = pd.Timestamp(now_dt)
    solpos = solarposition.get_solarposition(time=now, latitude=latitude, longitude=longitude)

    elevation = solpos['elevation'].iloc[0]
    azimuth = solpos['azimuth'].iloc[0]

    if elevation <= 0:
        return None

    print(f"  {now_dt.strftime('%I:%M %p')} - Elevation: {elevation:.2f}° Azimuth: {azimuth:.2f}°")

    azrad = np.deg2rad(azimuth)
    elrad = np.deg2rad(elevation)
    dx = np.sin(azrad)
    dy = -np.cos(azrad)
    tan_el = np.tan(elrad)

    with rasterio.open(TifPath) as src:
        resolution = src.res[0]
        dem = src.read(1)
        shadow = _cast_shadows(dem, src.height, src.width, dx, dy, tan_el, resolution, 25.0)
    return shadow

    
    
def create_water_mask(laz_file_path, ref_tif):
    """Creates a water mask by rasterizing class 9 (water) points onto the DEM grid."""
    laz = laspy.read(laz_file_path)
    water_pts = laz.classification == 9

    if water_pts.sum() == 0:
        with rasterio.open(ref_tif) as src:
            return np.zeros((src.height, src.width), dtype=bool)

    wx = laz.x[water_pts]
    wy = laz.y[water_pts]

    with rasterio.open(ref_tif) as src:
        transform = src.transform
        mask = np.zeros((src.height, src.width), dtype=bool)
        # map point coordinates to pixel indices
        cols = ((wx - transform.c) / transform.a).astype(int)
        rows = ((wy - transform.f) / transform.e).astype(int)
        # filter to valid bounds
        valid = (rows >= 0) & (rows < src.height) & (cols >= 0) & (cols < src.width)
        mask[rows[valid], cols[valid]] = True
    return mask

def render_shadow_frame(dem, shadow, structure_mask, water_mask, title):
    """Renders a frame with distinct colors: green=sunlit, dark blue=shadow, red=buildings/trees."""
    # normalize elevation to 0-1 for base color
    valid = dem[~np.isnan(dem)]
    vmin, vmax = np.min(valid), np.max(valid)
    norm = (dem - vmin) / (vmax - vmin + 1e-6)

    # build RGB image
    rgb = np.zeros((*dem.shape, 3))

    # sunlit ground: green tones
    rgb[:, :, 0] = 0.3 * norm          # slight red for warmth
    rgb[:, :, 1] = 0.5 + 0.4 * norm    # green
    rgb[:, :, 2] = 0.2 * norm          # slight blue

    # shadows: dark blue/purple
    if shadow is not None:
        shadow_mask = shadow == 1
        rgb[shadow_mask, 0] = 0.15
        rgb[shadow_mask, 1] = 0.1
        rgb[shadow_mask, 2] = 0.35

    # buildings/vegetation: red/orange
    if structure_mask is not None:
        rgb[structure_mask, 0] = 0.8
        rgb[structure_mask, 1] = 0.25
        rgb[structure_mask, 2] = 0.15

    # water: blue
    if water_mask is not None:
        rgb[water_mask, 0] = 0.1
        rgb[water_mask, 1] = 0.3
        rgb[water_mask, 2] = 0.7

    rgb = np.clip(rgb, 0, 1)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(rgb)
    ax.set_title(title, fontsize=14)
    ax.axis('off')

    # legend
    from matplotlib.patches import Patch
    legend_items = [
        Patch(facecolor=(0.5, 0.7, 0.3), label='Sunlit Ground'),
        Patch(facecolor=(0.15, 0.1, 0.35), label='Shadow'),
        Patch(facecolor=(0.8, 0.25, 0.15), label='Buildings / Trees'),
        Patch(facecolor=(0.1, 0.3, 0.7), label='Water'),
    ]
    ax.legend(handles=legend_items, loc='lower right', fontsize=10)
    return fig

def render_and_fill(laz_file, out_dir, name, exclude_cls):
    """Renders a DEM from LiDAR and fills holes. Returns path to clean DEM."""
    raw = os.path.join(out_dir, f"{name}_raw.tif")
    filled = os.path.join(out_dir, f"{name}_filled.tif")
    clean = os.path.join(out_dir, f"{name}_clean.tif")
    RenderLidar(laz_file, out_dir, f"{name}_raw.tif", exclude_cls, returns="all")
    fill_dem_holes(raw, filled)
    fill_dem_gaps(filled, clean)
    return clean

def generate_shadow_gif(shadow_tif, display_dem, structure_mask, water_mask, latitude, longitude, label, interpolated_directory):
    """Generates a shadow timelapse GIF using precomputed DEMs and masks."""
    cwd = os.getcwd()
    frames_dir = os.path.join(interpolated_directory, f"frames_{label}")
    os.makedirs(frames_dir, exist_ok=True)

    print(f"\n  [{label}] Casting shadows...")
    eastern = pytz.timezone('America/New_York')
    today = datetime.now(eastern).date()
    frames = []

    # calculate sunrise and sunset
    sun_times = solarposition.sun_rise_set_transit_ephem(
        pd.DatetimeIndex([pd.Timestamp(today, tz=eastern)]), latitude, longitude
    )
    sunrise = sun_times['sunrise'].iloc[0].astimezone(eastern)
    sunset = sun_times['sunset'].iloc[0].astimezone(eastern)
    print(f"  Sunrise: {sunrise.strftime('%I:%M %p')}  Sunset: {sunset.strftime('%I:%M %p')}")

    # start 1 hour after sunrise, end 1 hour before sunset
    start = eastern.localize(datetime(today.year, today.month, today.day, sunrise.hour + 1, 0, 0))
    end = eastern.localize(datetime(today.year, today.month, today.day, sunset.hour - 1, 0, 0))
    current = start

    while current <= end:
        time_dt = current
        shadow = ground_level_shadow_cast(shadow_tif, latitude, longitude, time_dt=time_dt)

        fig = render_shadow_frame(display_dem, shadow, structure_mask, water_mask,
                                  f"{label} - {time_dt.strftime('%I:%M %p')}")
        frame_path = os.path.join(frames_dir, f"frame_{time_dt.strftime('%H_%M')}.png")
        fig.savefig(frame_path, dpi=100, bbox_inches='tight')
        plt.close(fig)
        frames.append(frame_path)

        current = current + timedelta(minutes=15)

    # build gif in dated folder
    gif_dir = os.path.join(cwd, "output", "gifs", today.strftime("%Y-%m-%d"))
    os.makedirs(gif_dir, exist_ok=True)
    images = [imageio.imread(f) for f in frames]
    gif_path = os.path.join(gif_dir, f"{label}_timelapse.gif")
    imageio.mimsave(gif_path, images, duration=0.4, loop=0)
    print(f"  GIF saved to {gif_path}")

    # clean up frames
    shutil.rmtree(frames_dir)

def main():
    latitude = 33.96280751837907
    longitude = -83.38712753678465

    cwd = os.getcwd()
    interpolated_directory = os.path.join(cwd, "output")

    print("[1/6] Downloading LAZ file...")
    laz_file = os.path.join(cwd, fetch_lidar_file(latitude, longitude))

    print("[2/6] Clipping to radius...")
    laz_file = clip_laz_to_radius(laz_file, latitude, longitude, radius_ft=300)

    # render shared DEMs once
    print("[3/6] Rendering bare ground DEM...")
    bare_tif = render_and_fill(laz_file, interpolated_directory, "bare", "0,1,3,4,5,6")

    print("[4/6] Rendering ground + buildings DEM...")
    buildings_tif = render_and_fill(laz_file, interpolated_directory, "buildings", "0,1,3,4,5")

    print("[5/6] Rendering ground + buildings + trees DEM...")
    trees_tif = render_and_fill(laz_file, interpolated_directory, "trees", "0,1,3,4")

    # load shared arrays
    with rasterio.open(bare_tif) as src:
        bare_dem = src.read(1)
    with rasterio.open(buildings_tif) as src:
        display_dem = src.read(1)

    structure_mask = (display_dem - bare_dem) > 2.0
    water_mask = create_water_mask(laz_file, buildings_tif)

    # generate GIFs using shared DEMs
    print("[6/6] Generating GIFs...")
    generate_shadow_gif(buildings_tif, display_dem, structure_mask, water_mask,
                        latitude, longitude, "ground_only", interpolated_directory)
    generate_shadow_gif(trees_tif, display_dem, structure_mask, water_mask,
                        latitude, longitude, "with_vegetation", interpolated_directory)



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


    