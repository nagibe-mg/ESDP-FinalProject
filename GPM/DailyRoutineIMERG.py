### IMPORT MODULES
import ControlFlowIMERG as cfIMERG
import os
import datetime
from datetime import timedelta

### PARAMETER TO CHANGE ###
start = "2014-05-01 00:00:00"
end = "2014-05-02 00:00:00"
lat_min, lat_max = 14.0, 33.0  
lon_min, lon_max = -118.0, -86.0
region = [lat_min, lat_max, lon_min, lon_max]

# def daily_routine(start, end, region):
#     start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
#     end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
    
#     current = start_time
#     StoreZARR = "ZARR/IMERG_PRECIP.zarr"
    
#     # DAILY LOOP
#     while current <= end_time:
#         # Format current date as string
#         date_str = current.strftime('%Y%m%d')
        
#         # Format for download function (needs full datetime string)
#         current_str = current.strftime("%Y-%m-%d %H:%M:%S")
#         next_day = current + timedelta(days=1)
#         next_str = next_day.strftime("%Y-%m-%d %H:%M:%S")
        
#         tmp_file = f"temp_{date_str}.nc"  # Changed to .nc since you're using to_netcdf
        
#         """
#         Processing chain: Download -> Regrid -> Save
#         """
#         print(f"\nStarting pipeline for: {current_str}")
        
#         try: 
#             print(f"Downloading IMERG data for {date_str}...")
#             cfIMERG.download_IMERG(current_str, next_str, tmp_file)
            
#             print("Regridding to HEALPix (NSIDE 8 and 16).") 
#             regridded_data = cfIMERG.regrid_to_healpix(tmp_file, region)
            
#             print(f"Appending to Zarr store: {StoreZARR}") 
#             cfIMERG.save_to_zarr(regridded_data, StoreZARR, date_str)
            
#             # Clean up the temporary NetCDF file to save space
#             if os.path.exists(tmp_file):
#                 os.remove(tmp_file)
                
#             print(f"Successfully processed {date_str}")
            
#         except Exception as e:
#             print(f"Failed processing {date_str}: {e}")
#             if os.path.exists(tmp_file):
#                 os.remove(tmp_file)
        
#         # Move to next day
#         current += timedelta(days=1)
    
#     print("\nPipeline complete!")

def daily_routine(start, end):
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
    
    current = start_time
    StoreZARR = "ZARR/IMERG_PRECIP.zarr"  # Will create _nside8.zarr and _nside16.zarr
    
    while current <= end_time:
        date_str = current.strftime('%Y%m%d')
        current_str = current.strftime("%Y-%m-%d %H:%M:%S")
        next_day = current + timedelta(days=1)
        next_str = next_day.strftime("%Y-%m-%d %H:%M:%S")
        
        tmp_file = f"temp_{date_str}.nc"
        
        print(f"\nStarting pipeline for: {current_str}")
        
        try: 
            print(f"Downloading IMERG data for {date_str}...")
            cfIMERG.download_IMERG(current_str, next_str, tmp_file)
            
            print("Regridding to HEALPix (NSIDE 8 and 16)") 
            regridded_data = cfIMERG.regrid_to_healpix_fast(tmp_file)
            
            print(f"Saving to Zarr stores...") 
            cfIMERG.save_to_zarr(regridded_data, StoreZARR, date_str)
            
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
                
            print(f"Successfully processed {date_str}")
            
        except Exception as e:
            print(f"Failed processing {date_str}: {e}")
            import traceback
            traceback.print_exc()
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        
        current += timedelta(days=1)
    
    print("\nPipeline complete!")

## Run daily routine
daily_routine(start, end)