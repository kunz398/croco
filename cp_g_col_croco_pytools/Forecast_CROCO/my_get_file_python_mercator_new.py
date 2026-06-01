from croco_tools_params import *
import copernicusmarine
from datetime import datetime, timedelta
import re

def get_mercator(mercator_type, vars, geom, date, info, outname):
    """
    Generalized function to download Mercator data using the new Copernicus Marine API.
    Simplified signature: removed unnecessary parameters (pathmotu, url, sid, pid).
    
    Parameters:
    - mercator_type: 1 for Mercator 1/12°, 2 for UK Met Office 1/4°, 4 for independent variable datasets.
    - vars: String of variables, e.g., ' --variable thetao --variable so'.
    - geom: [min_lon, max_lon, min_lat, max_lat, min_depth, max_depth].
    - date: [start_date_str, end_date_str] in ISO format.
    - info: [username, password].
    - outname: Output filename.
    """
    
    # Authenticate
    username = info[0] if len(info) > 0 else None
    password = info[1] if len(info) > 1 else None
    if username and password:
        print(f"Authenticating with username: {username}")
        copernicusmarine.login(username=username, password=password)
        print("Authentication successful")
    else:
        print("Warning: No credentials provided. Ensure COPERNICUSMARINE_USERNAME and COPERNICUSMARINE_PASSWORD are set.")
    
    # Determine dataset_id based on mercator_type
    if mercator_type == 1:
        # Mercator 1/12° (from croco_tools_params: GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS / global-analysis-forecast-phy-001-024)
        dataset_id = "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"
    elif mercator_type == 2:
        # UK Met Office 1/4° (from croco_tools_params: GLOBAL_ANALYSISFORECAST_PHY_CPL_001_015-TDS / MetO-GLO-PHY-CPL-dm-TEM)
        dataset_id = "cmems_mod_glo_phy_coupled_anfc_0.25deg_P1D-m"
    elif mercator_type == 4:
        # Independent datasets for specific variables
        parsed_vars = re.findall(r'--variable\s+(\w+)', vars) if vars else []
        if 'uo' in parsed_vars or 'vo' in parsed_vars:
            dataset_id = "cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i"  # Currents
        elif 'so' in parsed_vars:
            dataset_id = "cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i"  # Salinity
        elif 'thetao' in parsed_vars:
            dataset_id = "cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i"  # Temperature
        elif 'ssh' in parsed_vars or 'sea_surface_height' in parsed_vars or 'zos' in parsed_vars:
            dataset_id = "cmems_mod_glo_phy_anfc_merged-sl_PT1H-i"  # Sea surface height
        else:
            raise ValueError(f"No matching dataset for variables {parsed_vars} in mercator_type=4. Supported: uo/vo, so, thetao, ssh/zos.")
    else:
        raise ValueError(f"Unsupported mercator_type: {mercator_type}. Use 1 for Mercator 1/12°, 2 for UK Met Office 1/4°, or 4 for independent variables.")
    
    # Parse variables from string to list
    if vars:
        variables = re.findall(r'--variable\s+(\w+)', vars)
    else:
        variables = None
    
    # Map variable names if necessary
    if variables:
        if dataset_id == "cmems_mod_glo_phy_anfc_merged-sl_PT1H-i":
            variables = ['sea_surface_height' if var == 'zos' else var for var in variables]
        # Add other mappings if needed
    
    # Parse dates
    start_datetime = datetime.fromisoformat(date[0].replace('Z', '+00:00'))
    end_datetime = datetime.fromisoformat(date[1].replace('Z', '+00:00'))
    
    # Parse geometry
    minimum_longitude = float(geom[0])
    maximum_longitude = float(geom[1])
    minimum_latitude = float(geom[2])
    maximum_latitude = float(geom[3])
    minimum_depth = float(geom[4]) if len(geom) > 4 and geom[4] != '0' else None
    maximum_depth = float(geom[5]) if len(geom) > 5 and geom[5] != '0' else None
    
    # Prepare subset parameters
    subset_kwargs = {
        "dataset_id": dataset_id,
        "variables": variables,
        "start_datetime": start_datetime,
        "end_datetime": end_datetime,
        "minimum_longitude": minimum_longitude,
        "maximum_longitude": maximum_longitude,
        "minimum_latitude": minimum_latitude,
        "maximum_latitude": maximum_latitude,
        "output_filename": outname,
    }
    
    # Add depth if specified
    if minimum_depth is not None:
        subset_kwargs["minimum_depth"] = minimum_depth
    if maximum_depth is not None:
        subset_kwargs["maximum_depth"] = maximum_depth
    
    # Download the data
    print(f"Starting download for {dataset_id}...")
    copernicusmarine.subset(**subset_kwargs)
    print(f"Download completed: {outname}")

# Example usage (uncomment to test)
# get_mercator(
#     mercator_type=1,
#     vars=' --variable thetao --variable so --variable uo --variable vo',
#     geom=['-173', '-165', '-23', '-16', '0', '6000'],
#     date=['2023-01-01T00:00:00.000Z', '2023-01-02T00:00:00.000Z'],
#     info=['your_username', 'your_password'],
#     outname='test_mercator.nc'
# )
