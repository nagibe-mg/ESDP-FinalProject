import datetime
from datetime import timedelta
import os
import gpm
import datetime
import xarray as xr
import healpy as hp
import numpy as np
import zarr

## functions

def download_IMERG(start, end):
    
    # SETUP
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")

    # Specify the product and product type
    product = "IMERG-FR"  # FR = final run
    
    product_type = "RS"  # "NRT"
    storage = "GES_DISC"  # or PPS, but I used the NASA Earthdata login information
    # Specify the version
    version = 7

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


def loadIMERG(start, end, filename):
    # SETUP
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
    
    # Specify the product and product type
    product = "IMERG-FR"  # FR = final run
    product_type = "RS"  # "NRT"
    
    # Specify the version
    version = 7
    
    # Load the dataset using the same parameters
    ds = gpm.open_dataset(
        product=product,
        product_type=product_type,
        version=version,
        start_time=start_time,
        end_time=end_time,
        variables = ['precipitation']
    )
    
    ds = (ds['precipitation'] * 0.5).sum(dim='time')

    # FORCE DATA LOADING - This is the key fix!
    print("Loading data into memory...")
    ds = ds.load()  # This actually downloads the data


    
    # Save with your custom filename
    print(f"Saving to {filename}...")
    ds.to_netcdf(filename, compute=True)
    print(f"Data saved to: {filename}")
    
    return ds


def regrid_to_healpix(source_file, date_str, nside_val=128):
    """
    Regrid lat-lon data to HEALPix and preserve the date
    
    Args:
        source_file: Path to NetCDF file
        date_str: Date string in format 'YYYY-MM-DD'
        nside_val: HEALPix NSIDE parameter
    
    Returns:
        xr.Dataset with HEALPix gridded data
    """
    ds = xr.open_dataset(source_file)
    
    # Convert date string to numpy datetime64
    date = np.datetime64(date_str)
    
    # Identify lat/lon dimensions
    lat_dim = 'lat' if 'lat' in ds.dims else 'latitude'
    lon_dim = 'lon' if 'lon' in ds.dims else 'longitude'
    
    lats = ds[lat_dim].values
    lons = ds[lon_dim].values
    
    # Create coordinate grids
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    # Convert to HEALPix angles
    theta = np.radians(90.0 - lat_grid).flatten()
    phi = np.radians(lon_grid).flatten()
    
    # Get HEALPix pixel indices
    pix_indices = hp.ang2pix(nside_val, theta, phi)
    npix = hp.nside2npix(nside_val)
    
    # Process all data variables
    data_vars = {}
    for var_name in ds.data_vars:
        if var_name in [lat_dim, lon_dim, 'crsWGS84']:
            continue
        
        data = ds[var_name].values
        
        # Initialize HEALPix map
        healpix_map = np.full(npix, hp.UNSEEN, dtype=np.float32)
        data_flat = data.flatten()
        
        # Average values into HEALPix pixels
        valid = ~np.isnan(data_flat)
        if np.any(valid):
            valid_pix = pix_indices[valid]
            valid_data = data_flat[valid]
            counts = np.bincount(valid_pix, minlength=npix)
            sums = np.bincount(valid_pix, weights=valid_data, minlength=npix)
            mask = counts > 0
            healpix_map[mask] = sums[mask] / counts[mask]
        
        # Create DataArray with time dimension
        data_vars[var_name] = xr.DataArray(
            healpix_map[np.newaxis, :],  # Add time dimension
            dims=['time', 'healpix_pixel'],
            coords={
                'time': [date],
                'healpix_pixel': np.arange(npix)
            }
        )
    
    # Create Dataset
    ds_healpix = xr.Dataset(data_vars)
    ds_healpix.attrs['nside'] = nside_val
    ds_healpix.attrs['date'] = date_str
    
    ds.close()
    return ds_healpix

def save_to_zarr(ds, store_path, nside):
    """
    Save HEALPix data to Zarr store
    
    Args:
        ds: xarray.Dataset with HEALPix data
        store_path: Base path for Zarr store
        nside: NSIDE parameter for naming
    """
    # Create NSIDE-specific path
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



def daily_routine(start, end, nside=128):
    """
    Process IMERG data day by day: Download -> Load -> Regrid -> Save to Zarr
    
    Args:
        start: Start date string 'YYYY-MM-DD HH:MM:SS'
        end: End date string 'YYYY-MM-DD HH:MM:SS'
        nside: HEALPix NSIDE parameter (default 128)
    """
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
    
    current = start_time
    store_zarr = "ZARR/IMERG_PRECIP.zarr"
    
    # Create ZARR directory if it doesn't exist
    os.makedirs("ZARR", exist_ok=True)
    
    # DAILY LOOP
    while current <= end_time:
        # Format current date
        date_str_file = current.strftime('%Y%m%d')
        date_str_iso = current.strftime('%Y-%m-%d')
        
        # Format for download function
        current_str = current.strftime("%Y-%m-%d %H:%M:%S")
        next_day = current + timedelta(days=1)
        next_str = next_day.strftime("%Y-%m-%d %H:%M:%S")
        
        tmp_file = f"temp_{date_str_file}.nc"
        
        print(f"\n{'='*60}")
        print(f"Processing: {date_str_iso}")
        print(f"{'='*60}")
        
        try:
            # Step 1: Download
            print(f"[1/4] Downloading IMERG data...")
            download_IMERG(current_str, next_str)
            
            # Step 2: Load and aggregate
            print(f"[2/4] Loading and aggregating data...")
            loadIMERG(current_str, next_str, tmp_file)
            
            # Step 3: Regrid to HEALPix
            print(f"[3/4] Regridding to HEALPix (NSIDE={nside})...")
            regridded_data = regrid_to_healpix(tmp_file, date_str_iso, nside_val=nside)
            
            # Step 4: Save to Zarr
            print(f"[4/4] Saving to Zarr store...")
            save_to_zarr(regridded_data, store_zarr, nside)
            
            # Clean up temporary file
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
                print(f"    Cleaned up temporary file: {tmp_file}")
            
            print(f"‚úì Successfully processed {date_str_iso}")
            
        except Exception as e:
            print(f"‚úó Failed processing {date_str_iso}: {e}")
            import traceback
            traceback.print_exc()
            
            # Clean up on failure
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        
        # Move to next day
        current += timedelta(days=1)
    
    print(f"\n{'='*60}")
    print("Pipeline complete!")
    print(f"{'='*60}")
