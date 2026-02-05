import xarray as xr
import healpy as hp
import numpy as np
import pandas as pd

def compute_precip_stats(
    store_path,
    nside=128,
    lat_min=14.0, lat_max=33.0,
    lon_min=-118.0, lon_max=-86.0
):

    # Load data
    ds = xr.open_zarr(store_path, consolidated=True)
    precip = ds["precipitation"]

    # HEALPix geometry
    npix = hp.nside2npix(nside)
    pixels = np.arange(npix)
    theta, phi = hp.pix2ang(nside, pixels)

    lats = 90.0 - np.degrees(theta)
    lons = np.degrees(phi)
    lons = np.where(lons > 180, lons - 360, lons)

    # Regional mask
    mask = (
        (lats >= lat_min) & (lats <= lat_max) &
        (lons >= lon_min) & (lons <= lon_max)
    )

    print(f"Pixels in region: {mask.sum()}")

    # -------------------------------------------------
    # DAILY STATISTICS (explicit but correct)
    # -------------------------------------------------
    daily_records = []

    for t in range(precip.sizes["time"]):
        data = precip.isel(time=t).values[mask]

        valid = ~np.isnan(data)

        daily_records.append({
            "date": pd.Timestamp(precip.time.values[t]).date(),
            "min_mm":  np.nanmin(data[valid]),
            "max_mm":  np.nanmax(data[valid]),
            "mean_mm": np.nanmean(data[valid]),
        })

    df_daily = pd.DataFrame(daily_records)
    df_daily.to_csv("daily_precip_stats.csv", index=False)

    return df_daily


store_path = "ZARR/IMERG_PRECIP_nside128.zarr"

compute_precip_stats(
    store_path,
    nside=128,
    lat_min=14.0, lat_max=33.0,
    lon_min=-118.0, lon_max=-86.0
)