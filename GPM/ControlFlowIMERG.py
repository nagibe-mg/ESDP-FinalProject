import gpm
from gpm.utils.geospatial import (
    get_circle_coordinates_around_point,
    get_country_extent,
    )
import datetime
from datetime import timedelta
import xarray as xr
import healpy as hp
import numpy as np
import os
import zarr

def download_IMERG(start, end, filename):
    
    # SETUP
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")

    # Specify the product and product type
    product = "IMERG-FR"  # FR = final run
    product_type = "RS"  # "NRT"
    storage = "GES_DISC"  # or PPS, but I used the NASA Earthdata login information
    # Specify the version
    version = 7

    #current = start_time

    # DAILY LOOP
    #while current <= end_time:

    # Download the data
    gpm.download(
        product=product,
        product_type=product_type,
        version=version,
        start_time=start_time,
        end_time=end_time,
        storage=storage,
        force_download=False,
        verbose=True,
        progress_bar=True,
        check_integrity=False,
    )
    # Load the dataset using the same parameters
    ds = gpm.open_dataset(
        product=product,
        product_type=product_type,
        version=version,
        start_time=start_time,
        end_time=end_time,
    )
    
    # Save with your custom filename
    ds.to_netcdf(filename)
    
    print(f"Data saved to: {filename}")
    #current += timedelta(days=1)


# regrid to healpix

def regrid_to_healpix_fast(source_file, nside_list=[8, 16]):
    """
    Fast HEALPix regridding using direct pixel assignment (like ERA5 code)
    """
    ds = xr.open_dataset(source_file)
    
    # Handle dimension names
    if 'valid_time' in ds.dims:
        ds = ds.rename({'valid_time': 'time'})
    
    lat_dim = 'lat' if 'lat' in ds.dims else 'latitude'
    lon_dim = 'lon' if 'lon' in ds.dims else 'longitude'
    
    # Get coordinate arrays
    lats = ds[lat_dim].values
    lons = ds[lon_dim].values
    
    # Create meshgrid
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    # Convert to theta (colatitude) and phi (longitude) in radians
    theta = np.radians(90.0 - lat_grid).flatten()
    phi = np.radians(lon_grid).flatten()
    
    regridded_data = {}
    
    for nside in nside_list:
        print(f"  ... Processing NSIDE={nside}")
        
        # Get HEALPix pixel indices for all grid points
        pix_indices = hp.ang2pix(nside, theta, phi)
        npix = hp.nside2npix(nside)
        
        # Process each variable
        data_vars = {}
        
        for var_name in ds.data_vars:
            # Skip coordinate variables
            if var_name in [lat_dim, lon_dim, 'crsWGS84']:
                continue
            
            print(f"    Processing variable: {var_name}")
            
            data = ds[var_name].values  # Shape: (time, lat, lon)
            
            # Initialize HEALPix map
            healpix_map = np.full((data.shape[0], npix), hp.UNSEEN, dtype=np.float32)
            
            # Flatten spatial dimensions
            data_flat = data.reshape(data.shape[0], -1)  # (time, lat*lon)
            
            # For each timestep, bin values into HEALPix pixels
            for t in range(data.shape[0]):
                valid = ~np.isnan(data_flat[t])
                
                if np.any(valid):
                    # Count how many grid points fall in each pixel
                    counts = np.bincount(pix_indices[valid], minlength=npix)
                    
                    # Sum values in each pixel
                    sums = np.bincount(
                        pix_indices[valid], 
                        weights=data_flat[t, valid], 
                        minlength=npix
                    )
                    
                    # Average (only for pixels with data)
                    mask = counts > 0
                    healpix_map[t, mask] = sums[mask] / counts[mask]
            
            # Create DataArray
            data_vars[var_name] = xr.DataArray(
                healpix_map,
                dims=['time', 'healpix_pixel'],
                coords={
                    'time': ds['time'],
                    'healpix_pixel': np.arange(npix)
                }
            )
        
        # Create Dataset
        ds_healpix = xr.Dataset(data_vars)
        ds_healpix.attrs['nside'] = nside
        ds_healpix.attrs.update(ds.attrs)  # Copy original attributes
        
        # Verify data
        if 'precipitation' in ds_healpix:
            precip_data = ds_healpix['precipitation'].values
            n_valid = np.sum(precip_data != hp.UNSEEN)
            print(f"    Valid pixels: {n_valid} / {precip_data.size}")
            valid_data = precip_data[precip_data != hp.UNSEEN]
            if len(valid_data) > 0:
                print(f"    Data range: {np.min(valid_data):.3f} to {np.max(valid_data):.3f}")
        
        regridded_data[nside] = ds_healpix
    
    ds.close()
    return regridded_data

# def regrid_to_healpix(source_file, nside_list=[8, 16]):
#     ds = xr.open_dataset(source_file)
    
#     if 'valid_time' in ds.dims:
#         ds = ds.rename({'valid_time': 'time'})
    
#     # Check dimension names
#     lat_dim = 'lat' if 'lat' in ds.dims else 'latitude'
#     lon_dim = 'lon' if 'lon' in ds.dims else 'longitude'
    
#     regridded_data = {}
    
#     for nside in nside_list:
#         print(f"  ... Processing NSIDE={nside}")
        
#         # Healpix grid
#         npix = hp.nside2npix(nside)
#         theta, phi = hp.pix2ang(nside, np.arange(npix))
        
#         # Convert to degrees and lat/long
#         target_lats = 90.0 - np.degrees(theta)
#         target_lons = np.degrees(phi)
        
#         # Structure as xarray DataArrays
#         target_lats_da = xr.DataArray(target_lats, dims="healpix_pixel")
#         target_lons_da = xr.DataArray(target_lons, dims="healpix_pixel")
        
#         # INTERPOLATE
#         ds_healpix = ds.interp(
#             {lat_dim: target_lats_da, lon_dim: target_lons_da},
#             method="nearest"  # Changed to nearest for speed
#         )
        
#         # *** CRITICAL: COMPUTE THE RESULT ***
#         ds_healpix = ds_healpix.compute()
        
#         # Verify we have data
#         if 'precipitation' in ds_healpix:
#             precip_data = ds_healpix['precipitation'].values
#             n_valid = np.sum(~np.isnan(precip_data))
#             print(f"    Valid pixels: {n_valid} / {precip_data.size}")
#             print(f"    Data range: {np.nanmin(precip_data):.3f} to {np.nanmax(precip_data):.3f}")
        
#         regridded_data[nside] = ds_healpix
    
#     ds.close()
#     return regridded_data

# saving to zarr store

# def save_to_zarr(data_dict, store_path, date_str):
#     for nside, ds in data_dict.items():
#         # Define Group Name (e.g., 'NSIDE_8')
#         group_name = f"NSIDE_{nside}"
        
#         # *** IMPORTANT: Ensure data is computed before saving ***
#         if hasattr(ds, 'compute'):
#             ds = ds.compute()
        
#         # Chunk by 1 day
#         chunks = {'time': 1, 'healpix_pixel': -1} 
#         ds = ds.chunk(chunks)
        
#         # Check if store exists
#         if not os.path.exists(store_path):
#             mode = 'w'
#             append_dim = None
#         else:
#             # Check if group exists
#             group_path = os.path.join(store_path, group_name)
#             if os.path.exists(group_path):
#                 mode = 'a'
#                 append_dim = 'time'
#             else:
#                 mode = 'a'  # Add new group to existing store
#                 append_dim = None
        
#         print(f"  Saving {group_name} to Zarr (mode={mode})...")
        
#         if append_dim:
#             ds.to_zarr(store_path, group=group_name, mode=mode, 
#                       append_dim=append_dim, consolidated=True)
#         else:
#             ds.to_zarr(store_path, group=group_name, mode=mode, 
#                       consolidated=True)
    
#     print(f"Data for {date_str} saved to {store_path}")

def save_to_zarr(data_dict, store_path, date_str):
    """
    Save HEALPix data to Zarr store (ERA5-style)
    """
    for nside, ds in data_dict.items():
        zarr_path = f"{store_path.replace('.zarr', '')}_nside{nside}.zarr"
        
        print(f"  Saving NSIDE={nside} to {zarr_path}")
        
        npix = hp.nside2npix(nside)
        
        # Check if zarr store exists
        if os.path.exists(zarr_path):
            # Append mode
            ds.to_zarr(zarr_path, mode='a', append_dim='time')
            # Reconsolidate metadata
            zarr.consolidate_metadata(zarr_path)
            print(f"    Appended to existing store")
        else:
            # Create new store with chunking
            encoding = {}
            for var in ds.data_vars:
                encoding[var] = {'chunks': (1, npix)}  # Chunk by time
            
            ds.to_zarr(zarr_path, mode='w', encoding=encoding, consolidated=True)
            print(f"    Created new store")

    # Delete the corrupted Zarr store
import shutil
if os.path.exists("ZARR/IMERG_PRECIP.zarr"):
    shutil.rmtree("ZARR/IMERG_PRECIP.zarr")
    print("Deleted old Zarr store")