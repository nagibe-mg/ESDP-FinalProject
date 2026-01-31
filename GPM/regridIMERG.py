import xarray as xr
import healpy as hp
import numpy as np


def regrid_to_healpix(source_file, nside_list=[8, 16]):
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