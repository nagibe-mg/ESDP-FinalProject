import downloadIMERG as do
import loadIMERG as lo
import regridIMERG as re
import saveZARR as sa

import datetime
from datetime import timedelta
import os

## trigger download of IMERG data
### PARAMETER TO CHANGE ###
start = "2014-05-01 00:00:00"
end = "2014-05-02 00:00:00"
#lat_min, lat_max = 14.0, 33.0  
#lon_min, lon_max = -118.0, -86.0
#region = [lat_min, lat_max, lon_min, lon_max]

def daily_routine(start, end): #, region=None):
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")
    
    current = start_time
    StoreZARR = "ZARR/IMERG_PRECIP.zarr"
    
    # DAILY LOOP
    while current <= end_time:
        # Format current date as string
        date_str = current.strftime('%Y%m%d')
        
        # Format for download function (needs full datetime string)
        current_str = current.strftime("%Y-%m-%d %H:%M:%S")
        next_day = current + timedelta(days=1)
        next_str = next_day.strftime("%Y-%m-%d %H:%M:%S")
        
        tmp_file = f"temp_{date_str}.nc"  # Changed to .nc since you're using to_netcdf
        
        """
        Processing chain: Download -> Regrid -> Save
        """
        print(f"\nStarting pipeline for: {current_str}")
        
        try: 
            print(f"Downloading IMERG data for {date_str}...")
            do.download_IMERG(current_str, next_str)
            
            print(f"Loading IMERG data for {date_str}...")
            lo.loadIMERG(current_str, next_str, tmp_file)
            
            

            print("Regridding to HEALPix (NSIDE 8 and 16).") 
            #re.regrid_to_healpix(tmp_file)
            regridded_data = re.regrid_to_healpix(tmp_file)
            
            print(f"Appending to Zarr store: {StoreZARR}") 
            sa.save_to_zarr(regridded_data, StoreZARR)
            
            # Clean up the temporary NetCDF file to save space
            #if os.path.exists(tmp_file):
            #    os.remove(tmp_file)
                
            print(f"Successfully processed {date_str}")
            
        except Exception as e:
            print(f"Failed processing {date_str}: {e}")
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        
        # Move to next day
        current += timedelta(days=1)
    
    print("\nPipeline complete!")

daily_routine(start, end)