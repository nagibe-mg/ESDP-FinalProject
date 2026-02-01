import os
import cdsapi
from datetime import datetime, timedelta
import numpy as np
import xarray as xr 
import healpy as hp
import matplotlib.pyplot as plt
from healpy.newvisufunc import projview
import json
import urllib.request

 
StoreZARR = "era5_precipitation.zarr"

# real routines

def check_if_processed(date):
    """
    checks if a specific date already exists in Zarr storage
    """
    return os.path.exists(f"processed_{date.strftime('%Y%m%d')}.txt")

def real_daily_data_routine(date,variable, stat):
    """
    Processing chain: Download -> Regrid -> Save
    """
    date_str = date.strftime('%Y-%m-%d')
    tmp_file = f"temp_{date_str}.nc"
    
    print(f"\nStarting pipeline for: {date.strftime('%Y-%m-%d')}")
    
    try: 
        
        print(f"Downloading {variable} for {date_str}...")
        download_era5(date, variable, stat, tmp_file)
        
        print("Regridding to HEALPix (NSIDE 128).") 
        regridded_data = regrid_to_healpix(tmp_file)
        
        
        print(f"Appending to Zarr store: {StoreZARR}") 
        save_to_zarr(regridded_data, StoreZARR, date_str)
    
        # mark as done
        with open(f"processed_{date.strftime('%Y%m%d')}.txt", "w") as f:
            f.write("done")
            
        #clean up the temporary NetCDF file to save space
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
            
        return True

    except Exception as e:
        print(f"Failed processing {date_str}: {e}")
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        return False

# download era5
def download_era5(date, variable, stat, targetFile):
    
    dataset = "derived-era5-single-levels-daily-statistics"
    request = {
        'product_type': 'reanalysis',
        'variable': variable,
        'year': date.strftime('%Y'),
        'month': date.strftime('%m'),
        'day': date.strftime('%d'),
        'daily_statistic': stat,
        'time_zone': 'utc+00:00',
        'frequency': '1_hourly',
        'area': [33, -118, 14, -86]
    }
    
    client = cdsapi.Client()
    client.retrieve(dataset, request, targetFile)


# regridding routine


def regrid_to_healpix(source_file, nside_list=[128]):
    ds = xr.open_dataset(source_file)
    
    if 'valid_time' in ds.dims:
        ds = ds.rename({'valid_time': 'time'})
    
    regridded_data = {}
    
    for nside in nside_list:
        print(f"  ... Processing NSIDE={nside}")
        
        # Helapix grid
        npix = hp.nside2npix(nside)
        # latitude and longitude for every pixel in radians: 
        #theta (colatitude) and phi (longitude) 
        theta, phi = hp.pix2ang(nside,np.arange(npix))
        
        # convert to degrees and lat/long
        target_lats = 90.0-np.degrees(theta)
        target_lons = np.degrees(phi)
        # align Longitude Conventions - convert our target (0-360) to match (-180 to 180).
        if ds.longitude.min() < 0:
            target_lons = np.where(target_lons > 180, target_lons - 360, target_lons)
        
        # INTERPOLATE
        # xarray's .interp(), pass the new target coordinates.
        #structure target_lats/lons as xarray DataArrays sharing a 'dim'
        target_lats_da = xr.DataArray(target_lats, dims="healpix_pixel")
        target_lons_da = xr.DataArray(target_lons, dims="healpix_pixel")
        
        ds_healpix = ds.interp(latitude=target_lats_da,longitude=target_lons_da,method="linear")
        
        regridded_data[nside] = ds_healpix

    ds.close()
    return regridded_data

# saving to zarr store

def save_to_zarr(data_dict, store_path, date_str):

    for nside, ds in data_dict.items():
        # Define Group Name (e.g., 'NSIDE_8')
        group_name = f"NSIDE_{nside}"
        # chunk by 1 day and keep spatial dims whole or split slightly.
        chunks = {'time': 1, 'healpix_pixel': -1} 
        ds = ds.chunk(chunks)
        
        #check if store exists, a mode for append and w for overwrite
        if not os.path.exists(store_path):
             mode = 'w'
             append_dim = None
        else:
             # check if group exist
             group_path = os.path.join(store_path, group_name)
             if os.path.exists(group_path):
                 mode = 'a'
                 append_dim = 'time'
             else:
                 mode = 'a' # Add new group to existing store
                 append_dim = None

        if append_dim:
            ds.to_zarr(store_path, group=group_name, mode=mode, append_dim=append_dim, consolidated=True)
        else:
            ds.to_zarr(store_path, group=group_name, mode=mode, consolidated=True)
            

    print(f"Data for {date_str} saved to {store_path}")
    
    
# plot routine


def add_mexico_border_standard():
    try:
        url = "https://raw.githubusercontent.com/johan/world.geo.json/master/countries/MEX.geo.json"
        with urllib.request.urlopen(url) as response:
            data = json.load(response)
        
        geom = data['features'][0]['geometry']
        coords_list = [geom['coordinates']] if geom['type'] == 'Polygon' else geom['coordinates']
        
        for poly in coords_list:
            ring = np.array(poly[0])
            plt.plot(ring[:, 0], ring[:, 1], 'k-', linewidth=1.5, alpha=1.0, zorder=10)
    except Exception as e:
        print(f"Could not load coastline: {e}")

def plot_precip_scatter(store_path, target_date, nside_val=128):
    nside_group = f"NSIDE_{nside_val}"
    ds_zarr = xr.open_zarr(store_path, group=nside_group, consolidated=True)
    
    try:
        sample = ds_zarr.sel(time=target_date, method='nearest')
    except KeyError:
        print("Date not found, using first index.")
        sample = ds_zarr.isel(time=0)

    # get Data
    precip_mm = sample['tp'].values * 1000
    precip_mm = np.nan_to_num(precip_mm, nan=0.0)
    print(f"Max Precip Value: {precip_mm.max():.2f} mm") 
    # coordinate calculation
    npix = len(precip_mm)
    theta, phi = hp.pix2ang(nside_val, np.arange(npix))
    lats = 90.0 - np.degrees(theta)
    lons = np.degrees(phi)
    lons = (lons + 180) % 360 - 180 # Force -180 to 180

    # Mexico region and min precip Box: [-118, -86, 14, 33]
    mask = (lons >= -120) & (lons <= -85) & (lats >= 12) & (lats <= 35) & (precip_mm > 0.05)
    
    plot_lons = lons[mask]
    plot_lats = lats[mask]
    plot_data = precip_mm[mask]

    # plot
    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(plot_lons, plot_lats, c=plot_data, cmap="Blues", 
                    vmin=0, vmax=60, s=65, marker='o', edgecolors='none')
    add_mexico_border_standard()
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.2, shrink=0.9)
    cbar.set_label("Precipitation (mm)", fontsize=12)

    ax.set_title(f"Total Precipitation | {target_date} | NSIDE {nside_val}", fontsize=14)
    ax.set_xlim([-118, -86])
    ax.set_ylim([14, 33])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.show()


# control flow

def real_run_control_flow(StartDate, EndDate, date_arg, variables, statistic):
    #date was provided
    if date_arg:
        desired_date = datetime.strptime(date_arg, "%Y-%m-%d")
        real_daily_data_routine(desired_date, variables, statistic)
        return

    # no arguments - look for oldest missing file 
    print(f"Scanning for missing data between {StartDate.date()} and {EndDate.date()}...")
    
    current_date = StartDate
    while current_date <= EndDate:
        if not check_if_processed(current_date):
            success = real_daily_data_routine(current_date, variables, statistic)
            
            if not success:
                print("Stopping --> Routine failed.")
                break 
        else:
            print(f"{current_date.date()} exists --> Skipping.")
        
        current_date += timedelta(days=1)