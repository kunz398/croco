from croco_tools_params import *
import copernicusmarine
import xarray as xr
from datetime import datetime, timedelta

# Define your parameters
dataset_id = "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i"  # Replace with the actual dataset ID
output_filename = "test_curr.nc"
variables = ["vo", "uo"] # Optional: specify variables to download
minimum_longitude = lonmin    # -173
maximum_longitude = lonmax    # -165
minimum_latitude = latmin     # -23
maximum_latitude = latmax     # -16

# Start datetime from croco_tools_param
start_datetime = datetime(Ymin, Mmin, Dmin, Hmin, Min_min, Smin)

# End datetime (start + hindcast + forecast days)  
end_datetime = start_datetime + timedelta(days=hdays + fdays)

# Download and subset the data
copernicusmarine.subset(
        dataset_id=dataset_id,
        output_filename=output_filename,
        variables=variables,
        start_datetime=start_datetime,
        end_datetime=end_datetime,
        minimum_latitude=minimum_latitude,
        maximum_latitude=maximum_latitude,
        minimum_longitude=minimum_longitude,
        maximum_longitude=maximum_longitude,
        # Add other relevant parameters as needed (e.g., depth, specific levels)
    )
