import gpm
import datetime

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
    )
    
    # FORCE DATA LOADING - This is the key fix!
    print("Loading data into memory...")
    ds = ds.load()  # This actually downloads the data
    
    # Save with your custom filename
    print(f"Saving to {filename}...")
    ds.to_netcdf(filename, compute=True)
    print(f"Data saved to: {filename}")
    
    return ds