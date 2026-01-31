import healpy as hp
import os
import zarr

def save_to_zarr(data_dict, store_path):
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