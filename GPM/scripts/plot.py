import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
import pandas as pd



def create_animation(store_path, start_idx=0, end_idx=None, nside=128,
                     lat_min=14.0, lat_max=33.0,
                     lon_min=-118.0, lon_max=-86.0,
                     output_file='animation.gif',
                     cmap='Blues', fps=1, dpi=300):
    """
    Create animation of daily precipitation
    
    Args:
        store_path: Path to Zarr store
        start_idx: Starting day index
        end_idx: Ending day index (None = all days)
        nside: HEALPix NSIDE parameter
        lat_min, lat_max, lon_min, lon_max: Region bounds
        output_file: Output filename (.gif or .mp4)
        cmap: Colormap
        fps: Frames per second
        dpi: Resolution
    """
    # Load data
    print("Loading data...")
    ds_zarr = xr.open_zarr(store_path, consolidated=True)
    npix = hp.nside2npix(nside)
    
    if end_idx is None:
        end_idx = len(ds_zarr.time)
    
    # Get pixel coordinates
    pixels = np.arange(npix)
    theta, phi = hp.pix2ang(nside, pixels)
    lats = 90.0 - np.degrees(theta)
    lons = np.where(np.degrees(phi) > 180, np.degrees(phi) - 360, np.degrees(phi))
    
    # Filter to region
    mask = (
        (lats >= lat_min) & (lats <= lat_max) &
        (lons >= lon_min) & (lons <= lon_max)
    )
    
    lats_region = lats[mask]
    lons_region = lons[mask]
    
    # Get data subset
    data_all = ds_zarr['precipitation'].isel(time=slice(start_idx, end_idx))
    num_frames = len(data_all.time)
    
    print(f"Creating animation for {num_frames} days...")
    
    # Calculate consistent color scale
    print("Calculating color scale...")
    all_values = []
    for i in range(num_frames):
        data_region = data_all.isel(time=i).values[mask]
        valid_mask = (data_region != hp.UNSEEN) & (~np.isnan(data_region))
        data_valid = data_region[valid_mask]
        all_values.extend(data_valid[data_valid > 0].tolist())
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10), 
                          subplot_kw={'projection': ccrs.PlateCarree()})
    
    cbar = None
    sc = None

    def update(frame):
        nonlocal sc, cbar

        ax.clear()

        # Get data for this frame
        daily_data = data_all.isel(time=frame)
        data_region = daily_data.values[mask]

        valid_mask = (data_region != hp.UNSEEN) & (~np.isnan(data_region))
        lats_plot = lats_region[valid_mask]
        lons_plot = lons_region[valid_mask]
        data_plot = data_region[valid_mask]

        sc = ax.scatter(
            lons_plot, lats_plot, c=data_plot,
            cmap=cmap, vmin=0, vmax=60,
            s=65, marker='o',
            transform=ccrs.PlateCarree(),
            edgecolors='none'
        )

        ax.coastlines(linewidth=1.5)
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        ax.set_extent([lon_min, lon_max, lat_min, lat_max])

        date_str = pd.Timestamp(daily_data.time.values).strftime('%Y-%m-%d')
        ax.set_title(
            f'Daily Precipitation | NSIDE={nside}\n{date_str}',
            fontsize=14, fontweight='bold'
        )

        # just one colorbar
        if cbar is None:
            cbar = fig.colorbar(
                sc, ax=ax,
                label='Precipitation (mm/day)',
                shrink=0.8, pad=0.05,
                orientation='horizontal'
            )

        print(f"  Frame {frame+1}/{num_frames}: {date_str}")
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=num_frames,
                        interval=1000/fps, repeat=True)
    
    # Save
    print(f"\nSaving animation to {output_file}...")
    if output_file.endswith('.gif'):
        writer = PillowWriter(fps=fps)
        anim.save(output_file, writer=writer, dpi=dpi)
    elif output_file.endswith('.mp4'):
        writer = FFMpegWriter(fps=fps, bitrate=1800)
        anim.save(output_file, writer=writer, dpi=dpi)
    else:
        raise ValueError("Output file must be .gif or .mp4")
    
    print(f"Animation saved: {output_file}")
    plt.close()
    
    return output_file


# ============================================================================
# USAGE EXAMPLES
# ============================================================================

if __name__ == "__main__":
    
    # Your region of interest
    lat_min, lat_max = 14.0, 33.0  
    lon_min, lon_max = -118.0, -86.0
    
    store_path = "../ZARR/IMERG_PRECIP_nside128.zarr"
    nside = 128
    
    print("="*60)
    print("IMERG Daily Precipitation Visualization")
    print("="*60)
    
    print("\n Creating animation...")
    create_animation(store_path, start_idx=0, end_idx=31, nside=nside,
                    lat_min=lat_min, lat_max=lat_max,
                    lon_min=lon_min, lon_max=lon_max,
                    output_file='../figures/precipitation_animation.gif',
                    fps=1, dpi=300)
    
    print("\n" + "="*60)
    print("All visualizations complete!")
    print("="*60)
