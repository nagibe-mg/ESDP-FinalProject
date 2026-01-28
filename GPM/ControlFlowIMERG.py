import gpm
from gpm.utils.geospatial import (
    get_circle_coordinates_around_point,
    get_country_extent,
    )
import datetime

def retrieve_and_process(start, end):
    
    # SETUP
    start_time = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end, "%Y-%m-%d %H:%M:%S")

    # Specify the product and product type
    product = "IMERG-FR"  # FR = final run
    product_type = "RS"  # "NRT"
    storage = "GES_DISC"  # or PPS, but I used the NASA Earthdata login information
    # Specify the version
    version = 7

    current = start_time

    # DAILY LOOP
    while current <= end_time:

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
    
    current += timedelta(days=1)