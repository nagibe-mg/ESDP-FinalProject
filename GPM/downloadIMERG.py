import gpm
import datetime

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
