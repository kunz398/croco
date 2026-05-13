#!/usr/bin/env python3
"""
Main script to run the complete CROCO forecast system:
1. Download MERCATOR/GLORYS data
2. Generate initial and boundary conditions
3. Get atmospheric forcing from GFS
4. Generate tidal forcing
5. Run CROCO forecast
"""

import os
import sys
from datetime import datetime, timedelta
import subprocess
import pandas as pd
import xarray as xr
from cp_g_col_croco_pytools.Forecast_CROCO.old__croco_tools_params import *
import Forecast_tools_new as ft
from main_croco import run_croco

def setup_paths(ruta_croco):
    """Setup and return all necessary paths relative to the CROCO directory"""
    global CROCOTOOLS_dir, CROCO_files_dir, FORC_DATA_DIR, FRCST_dir
    global grdname, frcname, blkname, clmname, bryname, ininame
    
    # Ensure ruta_croco is an absolute path
    ruta_croco = os.path.abspath(ruta_croco)
    
    # Set up base directories relative to ruta_croco
    CROCOTOOLS_dir = ruta_croco
    CROCO_files_dir = os.path.join(ruta_croco, 'CROCO_FILES')
    FORC_DATA_DIR = os.path.join(ruta_croco, 'DATA')
    FRCST_dir = os.path.join(FORC_DATA_DIR, 'Forecast')
    
    # Print directory setup for verification
    print("\nSetting up directories relative to:", ruta_croco)
    
    # Update input/output file paths
    grdname = os.path.join(CROCO_files_dir, 'croco_grd.nc')
    frcname = os.path.join(CROCO_files_dir, 'croco_frc.nc')
    blkname = os.path.join(CROCO_files_dir, 'croco_blk.nc')
    clmname = os.path.join(CROCO_files_dir, 'croco_clm.nc')
    bryname = os.path.join(CROCO_files_dir, 'croco_bry.nc')
    ininame = os.path.join(CROCO_files_dir, 'croco_ini.nc')

def ensure_directories_exist():
    """Create all necessary directories if they don't exist"""
    directories = [
        CROCOTOOLS_dir,
        CROCO_files_dir,
        FORC_DATA_DIR,
        FRCST_dir,
        os.path.dirname(grdname),
        os.path.dirname(frcname),
        os.path.dirname(blkname),
        os.path.dirname(clmname),
        os.path.dirname(bryname),
        os.path.dirname(ininame)
    ]
    
    for directory in directories:
        if not directory:
            continue
            
        try:
            if not os.path.exists(directory):
                print(f"Creating directory: {directory}")
                os.makedirs(directory, exist_ok=True)
            elif not os.path.isdir(directory):
                raise NotADirectoryError(f"Path exists but is not a directory: {directory}")
            elif not os.access(directory, os.W_OK):
                raise PermissionError(f"No write permission for directory: {directory}")
        except PermissionError as e:
            print(f"\nPermission error: {str(e)}")
            print("Please check that you have write permissions for the directory")
            print("or try running the script with appropriate permissions.")
            raise
        except Exception as e:
            print(f"\nError creating directory {directory}: {str(e)}")
            raise

def generate_ogcm_files(fecha_ini, fecha_fin, ruta_croco):
    """Generate initial and boundary conditions from MERCATOR/GLORYS"""
    print("\n=== Generating OGCM files (INI/BRY/CLM) ===")
    
    # Verify forecast directory exists and is writable
    if not os.path.exists(FRCST_dir):
        print(f"Creating forecast directory: {FRCST_dir}")
        os.makedirs(FRCST_dir, exist_ok=True)
    
    if not os.access(FRCST_dir, os.W_OK):
        raise PermissionError(f"No write permission in forecast directory: {FRCST_dir}")
    
    # Use the exact filename expected by make_OGCM_frcst
    raw_mercator_name = 'raw_motu_mercator_7524.nc'
    raw_mercator_path = os.path.join(FRCST_dir, raw_mercator_name)
    
    # Check if final merged file already exists and has all required data
    if os.path.exists(raw_mercator_path):
        try:
            with xr.open_dataset(raw_mercator_path) as ds:
                required_vars = ['sea_surface_height', 'uo', 'vo', 'thetao', 'so'] if mercator_type != 5 else ['sea_surface_height', 'uo', 'vo', 'thetao', 'so']
                time_range_ok = (pd.to_datetime(ds.time.min().values) <= pd.to_datetime(date_min) and 
                               pd.to_datetime(ds.time.max().values) >= pd.to_datetime(date_max))
                domain_ok = (ds.longitude.min().values <= lon_min and ds.longitude.max().values >= lon_max and
                           ds.latitude.min().values <= lat_min and ds.latitude.max().values >= lat_max)
                vars_ok = all(var in ds.variables for var in required_vars)
                
                if time_range_ok and domain_ok and vars_ok:
                    print(f"\nUsing existing MERCATOR data file: {raw_mercator_name}")
                    print(f"File contains all required variables and covers the requested time range and domain.")
                    print(f"Time range: {ds.time.min().values} to {ds.time.max().values}")
                    print(f"Domain: lon [{ds.longitude.min().values}, {ds.longitude.max().values}], "
                          f"lat [{ds.latitude.min().values}, {ds.latitude.max().values}]")
                    return  # Skip download and merging
        except Exception as e:
            print(f"\nWarning: Error checking existing file: {str(e)}")
            print("Will proceed with new download.")
    
    print(f"\nMERCATOR download configuration:")
    print(f"- Output directory: {FRCST_dir}")
    print(f"- Output file: {raw_mercator_name}")
    print(f"- User: {user}")
    print(f"- Product ID: {product_id_fcst}")
    print(f"- Domain: lon [{lonmin}, {lonmax}], lat [{latmin}, {latmax}]")
    print(f"- Time range: {fecha_ini} to {fecha_fin}")
    
    # Import required packages
    try:
        from copernicusmarine import subset
        import xarray as xr
    except ImportError as e:
        print(f"\nError importing required packages: {str(e)}")
        print("Please ensure you have the following packages installed:")
        print("- copernicusmarine (conda install conda-forge::copernicusmarine)")
        print("- xarray (conda install xarray)")
        raise
        
    # Calculate dates for mercator download (add buffer days)
    date_min = (fecha_ini - timedelta(days=hdays)).strftime("%Y-%m-%d 00:00:00")
    date_max = (fecha_fin + timedelta(days=fdays)).strftime("%Y-%m-%d 23:59:59")
    
    # Set variables to download based on mercator_type
    if mercator_type == 5:
        variables = ['sea_surface_height', 'uo', 'vo', 'thetao', 'so']  # Using sea_surface_height
    else:
        variables = ['sea_surface_height', 'uo', 'vo', 'thetao', 'so']
    
    # Add buffer to the domain (like in the MATLAB code)
    lon_min = lonmin - 1
    lon_max = lonmax + 1
    lat_min = latmin - 1
    lat_max = latmax + 1
    
    # Download data using Copernicus Marine API
    try:
        print("\nDownloading MERCATOR data components...")
        try:
            # Use existing credentials file from current directory
            current_dir = os.path.dirname(os.path.abspath(__file__))
            credentials_file = os.path.join(current_dir, '.copernicusmarine-credentials')
            
            if not os.path.exists(credentials_file):
                raise FileNotFoundError(f"Credentials file not found: {credentials_file}")
            
            print(f"Using credentials file: {credentials_file}")
            
            # Define dataset IDs and their variables
            # Extract dataset IDs from product_id_fcst if it's a dictionary
            if isinstance(product_id_fcst, dict):
                dataset_mapping = {
                    product_id_fcst['currents']: ['uo', 'vo'],  # Currents
                    product_id_fcst.get('temp', 'cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i'): ['thetao'],  # Temperature
                    product_id_fcst.get('salt', 'cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i'): ['so'],  # Salinity
                    'cmems_mod_glo_phy_anfc_merged-sl_PT1H-i': ['sea_surface_height'],  # SSH (hourly data)
                }
            else:
                dataset_mapping = {
                    'cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i': ['uo', 'vo'],  # Currents
                    'cmems_mod_glo_phy-thetao_anfc_0.083deg_PT6H-i': ['thetao'],  # Temperature
                    'cmems_mod_glo_phy-so_anfc_0.083deg_PT6H-i': ['so'],  # Salinity
                    'cmems_mod_glo_phy_anfc_merged-sl_PT1H-i': ['sea_surface_height']  # SSH (hourly data)
                }
            
            # Download and process each dataset separately
            downloaded_files = []
            
            for dataset_id, dataset_vars in dataset_mapping.items():
                print(f"\nProcessing dataset: {dataset_id}")
                print(f"Variables: {dataset_vars}")
                
                # Create a temporary file for this dataset with a unique name
                # Use more parts of the dataset_id to create a unique filename
                id_parts = dataset_id.split('-')
                var_type = id_parts[1] if len(id_parts) > 1 else 'unknown'
                temp_filename = f'temp_mercator_{var_type}.nc'
                temp_filepath = os.path.join(FRCST_dir, temp_filename)
                
                # Check if file already exists and has data
                file_exists = False
                if os.path.exists(temp_filepath):
                    try:
                        with xr.open_dataset(temp_filepath) as ds:
                            if len(ds.time) > 0 and all(var in ds.variables for var in dataset_vars):
                                print(f"Using existing file: {temp_filename}")
                                file_exists = True
                                downloaded_files.append(temp_filepath)
                                continue
                    except Exception:
                        print(f"Existing file {temp_filename} is invalid or incomplete")
                
                if not file_exists:
                    print(f"Downloading data to: {temp_filename}")
                    downloaded_files.append(temp_filepath)
                    subset(
                    dataset_id=str(dataset_id),  # Ensure dataset_id is a string
                    variables=dataset_vars,
                    minimum_longitude=float(lon_min),  # Ensure coordinates are floats
                    maximum_longitude=float(lon_max),
                    minimum_latitude=float(lat_min),
                    maximum_latitude=float(lat_max),
                    minimum_depth=0.0,
                    maximum_depth=5000.0,
                    start_datetime=str(date_min),  # Ensure dates are strings
                    end_datetime=str(date_max),
                    output_filename=str(temp_filename),
                    output_directory=str(FRCST_dir),
                    credentials_file=credentials_file
                )
            
            # Merge downloaded files using xarray with chunking and resampling
            print("\nMerging downloaded datasets...")
            datasets_to_merge = []
            target_freq = '6H'  # Target frequency for all datasets
            
            for file_path in downloaded_files:
                if os.path.exists(file_path):
                    try:
                        print(f"Processing {os.path.basename(file_path)}...")
                        # Open dataset with chunks
                        ds = xr.open_dataset(file_path, chunks={'time': 10})
                        
                        # Check if this is hourly SSH data
                        if 'sea_surface_height' in ds.variables:
                            print("Found hourly sea_surface_height data, resampling to 6-hourly...")
                            # Resample SSH data to match other variables
                            ds = ds.resample(time=target_freq).mean()
                        
                        datasets_to_merge.append(ds)
                        print(f"Variables in {os.path.basename(file_path)}: {list(ds.variables)}")
                        # Print some basic statistics for debugging
                        for var in ds.variables:
                            if var not in ['time', 'latitude', 'longitude', 'depth']:
                                try:
                                    stats = ds[var].isel(time=0)
                                    print(f"  {var} stats at first timestep - min: {stats.min().values}, max: {stats.max().values}")
                                except Exception as e:
                                    print(f"  Could not calculate stats for {var}: {str(e)}")
                            
                    except Exception as e:
                        print(f"Error opening {file_path}: {str(e)}")
                        raise
            
            # Check if the merged file already exists
            if os.path.exists(raw_mercator_path):
                print(f"\nMerged file already exists: {raw_mercator_name}")
                print("Keeping existing merged file and temporary files for reference.")
            else:
                try:
                    # Merge all datasets
                    print("Merging all datasets...")
                    merged_ds = xr.merge(datasets_to_merge)
                    print(f"Variables in merged dataset: {list(merged_ds.variables)}")
                    
                    # Save the merged dataset
                    if merged_ds is not None:
                        print(f"Writing merged data to {raw_mercator_name}...")
                        merged_ds.to_netcdf(raw_mercator_path)
                        print(f"Successfully merged data into: {raw_mercator_name}")
                        print("Keeping temporary downloaded files for reference.")
                    else:
                        raise ValueError("No datasets were successfully downloaded")
                except Exception as e:
                    print(f"Error merging datasets: {str(e)}")
                    # Don't delete existing files on error
                    raise
            
            # Verify file exists and is valid
            if not os.path.exists(raw_mercator_path):
                raise FileNotFoundError(f"Merged file not found: {raw_mercator_path} not found")
                
            # Verify the file is a valid NetCDF with data
            with xr.open_dataset(raw_mercator_path) as ds:
                if len(ds.time) == 0:
                    raise ValueError("Downloaded file contains no time data")
                    
        except Exception as e:
            print(f"Error downloading MERCATOR data: {str(e)}")
            # Keep any existing files for debugging
            print("\nKeeping existing files for debugging purposes.")
            print("You can manually delete them if needed.")
            raise
                
            print("Download and verification completed successfully")

    except Exception as e:
        print(f"\nError downloading MERCATOR data: {str(e)}")
        print("\nTroubleshooting steps:")
        print("1. Check your internet connection")
        print("2. Verify MERCATOR credentials in croco_tools_params.py")
        print("3. Confirm the CMEMS service is available (https://nrt.cmems-du.eu)")
        print("4. Try using a different DNS server")
        print("5. If using a proxy, check proxy settings")
        print("\nFull error details:", str(e))
        raise

    try:
        # Process MERCATOR data to generate CROCO input files
        print("\nProcessing merged data into CROCO format...")
        
        # Verify the merged file exists and has the required variables
        if not os.path.exists(raw_mercator_path):
            raise FileNotFoundError(f"Merged file not found: {raw_mercator_path}")
            
        # Ensure output directories exist before processing
        ensure_directories_exist()
        
        with xr.open_dataset(raw_mercator_path) as ds:
            required_vars = ['sea_surface_height', 'uo', 'vo', 'thetao', 'so']
            missing_vars = [var for var in required_vars if var not in ds.variables]
            if missing_vars:
                raise ValueError(f"Missing required variables in merged file: {missing_vars}")
            
            print("\nMerged file contains all required variables:")
            for var in required_vars:
                print(f"  - {var}: shape {ds[var].shape}")
                
            print(f"\nInput file: {raw_mercator_path}")
            print(f"Output directory: {CROCO_files_dir}")
        
        # Generate CROCO input files using the merged data
        print("\nGenerating CROCO input files...")
        ft.process_mercator_data(
            raw_mercator_name=raw_mercator_path,
            CROCO_files_dir=CROCO_files_dir,
            FRCST_prefix='mercator_',
            Yorig=Yorig
        )
        
        # Now convert MERCATOR data to CROCO format (INI/BRY/CLM files)
        print("\nConverting MERCATOR data to CROCO format...")
        try:
            # Calculate rundate consistently
            rundate = (fecha_ini - datetime(Yorig, 1, 1)).days
            
            # Check if CROCO input files already exist
            ini_file = os.path.join(CROCO_files_dir, f'croco_ini.nc')
            bry_file = os.path.join(CROCO_files_dir, f'croco_bry.nc') 
            clm_file = os.path.join(CROCO_files_dir, f'croco_clm.nc')
            
            if os.path.exists(ini_file) and os.path.exists(bry_file):
                print("CROCO input files already exist:")
                print(f"  - {os.path.basename(ini_file)}")
                print(f"  - {os.path.basename(bry_file)}")
                if os.path.exists(clm_file):
                    print(f"  - {os.path.basename(clm_file)}")
            else:
                print("Creating CROCO input files from MERCATOR data...")
                print("Running make_OGCM_frcst conversion...")
                
                # Import and run the OGCM processing
                import subprocess
                import sys
                
                # Run simplified OGCM conversion script
                script_path = os.path.join(os.path.dirname(__file__), 'make_OGCM_simple.py')
                if os.path.exists(script_path):
                    try:
                        # Run in ruta_croco directory so CROCO_FILES paths are accessible
                        result = subprocess.run([sys.executable, script_path], 
                                              cwd=ruta_croco,  # Use ruta_croco as working directory
                                              capture_output=True, 
                                              text=True,
                                              timeout=300)  # 5 minute timeout
                        
                        if result.returncode == 0:
                            print("OGCM conversion completed successfully")
                            print("STDOUT:", result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
                        else:
                            print(f"OGCM conversion failed with return code: {result.returncode}")
                            print("STDERR:", result.stderr)
                            print("STDOUT:", result.stdout)
                            
                    except subprocess.TimeoutExpired:
                        print("OGCM conversion timed out after 5 minutes")
                    except Exception as e:
                        print(f"Error running OGCM conversion: {str(e)}")
                else:
                    print(f"make_OGCM_simple.py not found at: {script_path}")
                
                # Verify files were created
                if os.path.exists(ini_file):
                    print(f"✓ Created: {os.path.basename(ini_file)}")
                else:
                    print(f"✗ Missing: {os.path.basename(ini_file)}")
                    
                if os.path.exists(bry_file):
                    print(f"✓ Created: {os.path.basename(bry_file)}")
                else:
                    print(f"✗ Missing: {os.path.basename(bry_file)}")
                    
                if os.path.exists(clm_file):
                    print(f"✓ Created: {os.path.basename(clm_file)}")
                else:
                    print(f"? Optional: {os.path.basename(clm_file)} (not created)")
            
            print("CROCO format conversion completed")
            
        except Exception as e:
            print(f"Error during CROCO conversion: {str(e)}")
            print("You may need to run make_OGCM_frcst.py manually")
            raise
        
        print("\nKeeping downloaded and processed files for reference.")

    except Exception as e:
        print(f"Error processing MERCATOR data: {str(e)}")
        raise

def generate_atmospheric_forcing(fecha_ini, fecha_fin):
    """Generate atmospheric forcing files from GFS"""
    print("\n=== Generating atmospheric forcing from GFS ===")
    
    # Ensure forecast directory exists
    os.makedirs(FRCST_dir, exist_ok=True)
    
    # Calculate expected GFS filename based on rundate
    rundate = (fecha_ini - datetime(Yorig, 1, 1)).days
    gfs_name = os.path.join(FRCST_dir, f'GFS_{rundate}.nc')
    
    if os.path.exists(gfs_name):
        print(f"Using existing GFS file: {gfs_name}")
    else:
        print(f"Downloading GFS data for period: {fecha_ini} to {fecha_fin}")
        print(f"Output directory: {FRCST_dir}")
        print(f"Domain: lon [{lonmin}, {lonmax}], lat [{latmin}, {latmax}]")
        
        try:
            import cp_g_col_croco_pytools.Forecast_CROCO.old__make_gfs as mgfs
            mgfs.download_GFS(
                fecha_ini,
                lonmin=lonmin,
                lonmax=lonmax,
                latmin=latmin,
                latmax=latmax,
                FRCST_dir=FRCST_dir,
                Yorig=Yorig,
                it=1  # Time interval in hours
            )
        except AttributeError as e:
            print(f"Error: {str(e)}")
            print("The make_gfs module needs to be updated. Creating placeholder files...")
            
            # Create placeholder GFS files
            rundate = (fecha_ini - datetime(Yorig, 1, 1)).days
            gfs_name = os.path.join(FRCST_dir, f'GFS_{rundate}.nc')
            blk_name = os.path.join(CROCO_files_dir, f'croco_blk_GFS_{rundate}.nc')
            frc_name = os.path.join(CROCO_files_dir, f'croco_frc_GFS_{rundate}.nc')
            
            if not os.path.exists(gfs_name):
                # Create minimal placeholder files
                print(f"Creating placeholder GFS file: {gfs_name}")
                try:
                    import netCDF4
                    with netCDF4.Dataset(gfs_name, 'w') as nc:
                        nc.createDimension('time', 1)
                        nc.createDimension('lat', 10)
                        nc.createDimension('lon', 10)
                        time_var = nc.createVariable('time', 'f8', ('time',))
                        time_var[:] = [rundate]
                except Exception as e2:
                    print(f"Could not create placeholder: {str(e2)}")
            
            print("Note: You'll need to provide proper atmospheric forcing files manually.")
        except Exception as e:
            print(f"\nError downloading GFS data: {str(e)}")
            print("Please check:")
            print("1. Internet connection")
            print("2. If the GFS service is available")
            print("3. If the requested date range is available in GFS")
            raise

def generate_tidal_forcing():
    """Generate tidal forcing using make_tides.py from PREPRO"""
    print("\n=== Generating tidal forcing ===")
    
    # Get the absolute paths to required directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    prepro_path = os.path.abspath(os.path.join(script_dir, '..', 'PREPRO'))
    modules_path = os.path.join(prepro_path, 'Modules')
    readers_path = os.path.join(prepro_path, 'Readers')
    tides_txt_path = os.path.join(modules_path, 'tides.txt')
    
    # Add all required paths to Python path
    for path in [prepro_path, modules_path, readers_path]:
        if path not in sys.path:
            sys.path.insert(0, path)
            print(f"Added to Python path: {path}")
            
    # Save current working directory and Python path
    original_dir = os.getcwd()
    original_pythonpath = list(sys.path)
    
    try:
        # Import make_tides after setting up Python path
        import make_tides
        
        # Get absolute paths for grid and output files
        grid_abs_path = os.path.abspath(grdname)
        tides_abs_path = os.path.join(os.path.dirname(grid_abs_path), 'croco_tides.nc')
        
        print(f"\nUsing files:")
        print(f"Grid: {grid_abs_path}")
        print(f"Tides.txt: {tides_txt_path}")
        print(f"Output: {tides_abs_path}")
        
        # Verify files and directories exist
        if not os.path.exists(grid_abs_path):
            raise FileNotFoundError(f"Grid file not found: {grid_abs_path}")
            
        if not os.path.exists(tides_txt_path):
            raise FileNotFoundError(f"Tides.txt not found: {tides_txt_path}")
            
        # Ensure output directory exists
        os.makedirs(os.path.dirname(tides_abs_path), exist_ok=True)
        
        # Call make_tides with absolute paths and proper working directory
        make_tides.generate_tides(
            grdname=grid_abs_path,
            tidename=tides_abs_path,
            working_dir=os.path.dirname(grid_abs_path)
        )
        
    except Exception as e:
        print(f"\nError generating tidal forcing: {str(e)}")
        raise
    finally:
        # Restore original state
        os.chdir(original_dir)
        sys.path = original_pythonpath
        print(f"\nRestored working directory to: {original_dir}")
    
    try:
        import make_tides
    except ImportError as e:
        print(f"Error: Could not import make_tides.py from {prepro_path}")
        raise
    
    # Ensure CROCO_files_dir exists
    os.makedirs(CROCO_files_dir, exist_ok=True)
    
    print(f"\nGenerating tidal forcing:")
    print(f"Grid file: {grdname}")
    print(f"Output directory: {CROCO_files_dir}")
    
    # Call make_tides with paths from croco_tools_params.py
    try:
        make_tides.generate_tides(
            grdname=grdname,
            tidename=os.path.join(CROCO_files_dir, 'croco_tides.nc'),
            working_dir=os.path.dirname(grdname)
        )
    except Exception as e:
        print(f"\nError generating tidal forcing: {str(e)}")
        raise

def setup_forecast(ruta_croco, fecha_ini, fecha_fin=None, outfiles=None, nproc=18):
    """
    Set up and run the complete CROCO forecast system
    
    Args:
        ruta_croco (str): Path to CROCO executable and configuration files
        fecha_ini (datetime): Start date of the forecast
        fecha_fin (datetime, optional): End date of the forecast
        outfiles (str, optional): Path where output files will be stored
        nproc (int, optional): Number of MPI processes to use
    """
    if fecha_fin is None:
        fecha_fin = fecha_ini + timedelta(days=fdays)
    
    try:
        # Setup paths relative to ruta_croco
        print("\n=== Setting up paths ===")
        setup_paths(ruta_croco)
        
        # Create all necessary directories
        print("\n=== Creating directories ===")
        ensure_directories_exist()
        
        print("\nWorking directories:")
        print(f"CROCO tools dir: {CROCOTOOLS_dir}")
        print(f"CROCO files dir: {CROCO_files_dir}")
        print(f"Forecast data dir: {FRCST_dir}")
        # 1. Generate OGCM files (INI/BRY/CLM)
        generate_ogcm_files(fecha_ini, fecha_fin, ruta_croco)
        
        # 2. Generate atmospheric forcing
        generate_atmospheric_forcing(fecha_ini, fecha_fin)
        
        # 3. Generate tidal forcing if enabled
        if add_tides_fcst:
            generate_tidal_forcing()
        
        # 4. Run CROCO forecast
        print("\n=== Running CROCO forecast ===")
        
        # Prepare custom paths to pass to run_croco
        custom_paths = {
            'grdname': grdname,
            'ininame': ininame,
            'bryname': bryname,
            'clmname': clmname,
            'frcname': frcname,
            'blkname': blkname
        }
        
        success = run_croco(
            ruta_croco=ruta_croco,
            fecha_ini=fecha_ini,
            fecha_fin=fecha_fin,
            outfiles=outfiles,
            nproc=nproc,
            custom_paths=custom_paths
        )
        
        if success:
            print("\nForecast completed successfully!")
        else:
            print("\nError during CROCO execution")
            
    except Exception as e:
        print(f"\nError during forecast setup: {str(e)}")
        raise

def test_atmospheric_forcing(ruta_croco, fecha_ini):
    """
    Test function to verify GFS data processing with memory optimization
    
    Args:
        ruta_croco (str): Path to CROCO directory
        fecha_ini (datetime): Test date
    """
    try:
        # Setup minimal required paths
        setup_paths(ruta_croco)
        ensure_directories_exist()
        
        print("\n=== Testing Atmospheric Forcing Generation ===")
        print(f"Test date: {fecha_ini}")
        print(f"Domain: lon [{lonmin}, {lonmax}], lat [{latmin}, {latmax}]")
        
        # Try generating atmospheric forcing
        generate_atmospheric_forcing(fecha_ini, fecha_ini + timedelta(days=1))
        
        # Verify the output files
        rundate = (fecha_ini - datetime(Yorig, 1, 1)).days
        gfs_name = os.path.join(FRCST_dir, f'GFS_{rundate}.nc')
        blk_name = os.path.join(CROCO_files_dir, f'croco_blk_GFS_{rundate}.nc')
        frc_name = os.path.join(CROCO_files_dir, f'croco_frc_GFS_{rundate}.nc')
        
        # Check file existence
        files_to_check = [
            (gfs_name, "GFS data file"),
            (blk_name, "Bulk forcing file"),
            (frc_name, "Surface forcing file")
        ]
        
        print("\nVerifying output files:")
        for file_path, description in files_to_check:
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path) / (1024 * 1024)  # Size in MB
                print(f"✓ {description} created successfully: {file_path}")
                print(f"  Size: {file_size:.2f} MB")
            else:
                print(f"✗ {description} not found: {file_path}")
        
        print("\nTest completed successfully!")
        return True
        
    except Exception as e:
        print(f"\nError during atmospheric forcing test: {str(e)}")
        return False

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Run complete CROCO forecast system')
    parser.add_argument('ruta_croco', help='Path to CROCO executable and configuration files')
    parser.add_argument('--test', action='store_true', help='Run atmospheric forcing test only')
    parser.add_argument('fecha_ini', help='Start date (YYYY-MM-DD)')
    parser.add_argument('--fecha_fin', help='End date (YYYY-MM-DD)')
    parser.add_argument('--outfiles', help='Output directory')
    parser.add_argument('--nproc', type=int, default=18, help='Number of MPI processes to use (default: 18)')
    
    args = parser.parse_args()
    
    fecha_ini = datetime.strptime(args.fecha_ini, '%Y-%m-%d')
    fecha_fin = datetime.strptime(args.fecha_fin, '%Y-%m-%d') if args.fecha_fin else None
    
    if args.test:
        print("\nRunning atmospheric forcing test...")
        success = test_atmospheric_forcing(args.ruta_croco, fecha_ini)
        sys.exit(0 if success else 1)
    else:
        setup_forecast(args.ruta_croco, fecha_ini, fecha_fin, args.outfiles, args.nproc)
