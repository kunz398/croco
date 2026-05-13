#!/usr/bin/env python3
"""
CROCO Forecast System - Python Version
=====================================

This script translates the original MATLAB-based run_croco_forecast.bash 
to work with the Python forecast tools. It maintains the same functionality
and workflow while leveraging the Python implementations.

Author: Translated from MATLAB version
Date: August 2025
"""

from html import parser
import os
import sys
import shutil
import subprocess
import datetime
import argparse
from pathlib import Path
import logging

# Import CROCO Python tools
try:
    from croco_tools_params import *
    # Don't import make_gfs and make_OGCM_frcst directly to avoid side effects
    # import make_gfs
    # import make_OGCM_frcst
    import Preprocessing_tools as ppt
    import Forecast_tools as ft
except ImportError as e:
    print(f"Error importing CROCO tools: {e}")
    print("Make sure all Python CROCO tools are in the same directory")
    sys.exit(1)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CROCOForecastRunner:
    """
    Main class for running CROCO forecasts using Python tools
    """
    
    def __init__(self, config_file=None):
        """Initialize the forecast runner with configuration"""
        
        # User-defined parameters (equivalent to bash script parameters)
        self.CLEAN = 0  # Clean results of previous forecasts
        self.PRE_PROCESS = 1  # Get forcing files and process them
        self.RESTART = 0  # Restart from previous forecast
        self.RUN = 1  # Run hindcast/forecast
        self.PLOT = 0  # Make plots
        
        # Time parameters
        self.DT = 160  # Model time step [seconds]
        self.NFAST = 90  # Number of barotropic time steps within one baroclinic
        self.NDAYS_HIND = 1  # Hindcast depth [days]
        self.NDAYS_FCST = 7  # Forecast depth [days]
        self.ND_AVG = 24  # Output frequency [hours] - average
        self.ND_HIS = 1   # Output frequency [hours] - history
        self.ND_RST = 24  # Output frequency [hours] - restart
        
        # Execution parameters (Linux-focused) - Define BEFORE using in setup_file_names
        self.EXEC = "mpirun -np 18"  # Adjust based on your system
        self.MODEL = "croco"
        self.CODFILE = "./croco"
        
        self.datestr = datetime.datetime.now().strftime('%Y-%m-%d')
        
        # Directories and files - setup AFTER defining MODEL
        self.setup_directories()
        self.setup_file_names()
        
        # Load tidal forcing configuration
        try:
            from croco_tools_params import add_tides_fcst, grdname
            self.add_tides_fcst = add_tides_fcst
            self.grdname = grdname
            logger.info(f"Tidal forcing configuration: add_tides_fcst = {self.add_tides_fcst}")
            logger.info(f"Grid file: {self.grdname}")
        except ImportError:
            self.add_tides_fcst = 0
            self.grdname = None
            logger.warning("add_tides_fcst or grdname not found in croco_tools_params, defaulting to disabled")
        
    def setup_directories(self):
        """Setup directory structure"""
        # Use RUN_dir from croco_tools_params if available, otherwise use current directory
        try:
            from croco_tools_params import RUN_dir
            self.RUNDIR = RUN_dir
        except ImportError:
            self.RUNDIR = os.getcwd()
            
        self.INPUTDIR = self.RUNDIR
        self.SCRATCHDIR = os.path.join(self.RUNDIR, "SCRATCH")
        self.MSSDIR = os.path.join(self.RUNDIR, "CROCO_FILES")
        self.MSSOUT = os.path.join(self.RUNDIR, "FORECAST")
        
        # Create directories if they don't exist
        for directory in [self.SCRATCHDIR, self.MSSOUT]:
            os.makedirs(directory, exist_ok=True)
            
        logger.info(f"Using RUN directory: {self.RUNDIR}")
        logger.info(f"CROCO_FILES directory: {self.MSSDIR}")
        logger.info(f"SCRATCH directory: {self.SCRATCHDIR}")
        logger.info(f"FORECAST output directory: {self.MSSOUT}")
            
    def setup_file_names(self):
        """Setup file names following CROCO conventions"""
        self.GRDFILE = f"{self.MODEL}_grd.nc"
        self.INIFILE = f"{self.MODEL}_ini.nc"
        self.RSTFILE = f"{self.MODEL}_rst.nc"
        self.AVGFILE = f"{self.MODEL}_avg.nc"
        self.HISFILE = f"{self.MODEL}_his.nc"
        self.BLKFILE = f"{self.MODEL}_blk_GFS_0.nc"
        self.FRCFILE = f"{self.MODEL}_frc_GFS_0.nc"
        self.BRYFILE = f"{self.MODEL}_bry_mercator_0.nc"
        self.CLMFILE = f"{self.MODEL}_clm_mercator_0.nc"
        self.INIFILE0 = f"{self.MODEL}_ini_mercator_0.nc"
        
    def clean_previous_results(self):
        """Clean results from previous forecasts"""
        if self.CLEAN:
            logger.info("Cleaning previous forecast results...")
            
            # Clean forcing files
            patterns = [
                f"{self.MODEL}_blk_*",
                f"{self.MODEL}_frc_*", 
                f"{self.MODEL}_bry_*",
                f"{self.MODEL}_clm_*",
                f"{self.MODEL}_ini_*"
            ]
            
            for pattern in patterns:
                for file in Path(self.MSSDIR).glob(pattern):
                    file.unlink()
                    
            # Clean output files
            output_patterns = [
                f"{self.MODEL}_his_*",
                f"{self.MODEL}_avg_*",
                f"{self.MODEL}_rst_*"
            ]
            
            for pattern in output_patterns:
                for file in Path(self.MSSOUT).glob(pattern):
                    file.unlink()
                    
            logger.info("Previous results cleaned")
            
    def preprocess_forcing(self):
        """Process boundary and forcing files using Python tools, with correct skip/only GFS logic"""
        if not self.PRE_PROCESS:
            logger.info("Skipping preprocessing (PRE_PROCESS=0)")
            return True
        logger.info("Processing boundary and forcing files...")
        try:
            # Check if parameter file exists
            if not os.path.exists('croco_tools_params.py'):
                logger.error("croco_tools_params.py not found. Please ensure it exists.")
                return False
            # Test if we can import the parameters
            try:
                import croco_tools_params
                logger.info("Configuration loaded successfully")
                self.makeblk_gfs = getattr(croco_tools_params, 'makeblk_gfs', 1)
                self.makefrc_gfs = getattr(croco_tools_params, 'makefrc_gfs', 1)
                logger.info(
                    f"GFS atmospheric output switches: makeblk_gfs={self.makeblk_gfs}, makefrc_gfs={self.makefrc_gfs}"
                )
            except Exception as e:
                logger.error(f"Error loading configuration: {e}")
                return False
            # Determine which processing steps to run
            if getattr(self, 'SKIP_GFS', False):
                run_gfs = False
                run_ogcm = True
            elif getattr(self, 'ONLY_GFS', False):
                run_gfs = True
                run_ogcm = False
            else:
                run_gfs = True
                run_ogcm = True
            # If both skip_gfs and only_gfs are set, skip all (user error)
            if getattr(self, 'SKIP_GFS', False) and getattr(self, 'ONLY_GFS', False):
                logger.error("Both --skip-gfs and --only-gfs flags set. Nothing to process.")
                return False
            # Run GFS processing if requested
            if run_gfs:
                logger.info("Running GFS data processing...")
                try:
                    result = subprocess.run([sys.executable, 'make_gfs.py'], 
                                            capture_output=True, text=True)
                    if result.returncode != 0:
                        logger.error(f"GFS processing failed: {result.stderr}")
                        if result.stdout:
                            logger.error(f"GFS stdout: {result.stdout}")
                        return False
                    logger.info("GFS processing completed")
                except Exception as e:
                    logger.error(f"Error running GFS processing: {e}")
                    return False
            else:
                logger.info("Skipping GFS data processing (--skip-gfs)")
            # Run OGCM/Mercator processing if requested
            if run_ogcm:
                logger.info("Running OGCM (Mercator) data processing...")
                try:
                    result = subprocess.run([sys.executable, 'make_OGCM_frcst.py'],
                                            capture_output=True, text=True)
                    if result.returncode != 0:
                        logger.error(f"OGCM processing failed: {result.stderr}")
                        if result.stdout:
                            logger.error(f"OGCM stdout: {result.stdout}")
                        return False
                    logger.info("OGCM processing completed")
                except Exception as e:
                    logger.error(f"Error running OGCM processing: {e}")
                    return False
            else:
                logger.info("Skipping OGCM (Mercator) data processing (--only-gfs)")
            return True
        except Exception as e:
            logger.error(f"Unexpected error during preprocessing: {e}")
            return False
            
    def copy_files_to_scratch(self):
        """Copy necessary files to scratch directory"""
        logger.info("Copying files to scratch directory...")
        
        # File patterns to copy from MSSDIR (use glob patterns to find most recent files)
        makeblk_gfs = getattr(self, 'makeblk_gfs', 1)
        makefrc_gfs = getattr(self, 'makefrc_gfs', 1)

        file_patterns = [
            (f"{self.MODEL}_bry_mercator_*.nc", self.BRYFILE),
            (f"{self.MODEL}_clm_mercator_*.nc", self.CLMFILE),
            (f"{self.MODEL}_grd.nc", self.GRDFILE)
        ]

        if makeblk_gfs:
            file_patterns.insert(0, (f"{self.MODEL}_blk_GFS_*.nc", self.BLKFILE))
        else:
            logger.info("Skipping BLK copy (makeblk_gfs=0)")

        if makefrc_gfs:
            insert_at = 1 if makeblk_gfs else 0
            file_patterns.insert(insert_at, (f"{self.MODEL}_frc_GFS_*.nc", self.FRCFILE))
        else:
            logger.info("Skipping FRC copy (makefrc_gfs=0)")
        
        # Add tidal forcing file if tidal forcing is enabled
        if self.add_tides_fcst == 1:
            tidal_pattern = f"{self.MODEL}_frc.nc"
            file_patterns.append((tidal_pattern, f"{self.MODEL}_frc.nc"))
        
        for pattern, dst_filename in file_patterns:
            # Find the most recent file matching the pattern
            matching_files = list(Path(self.MSSDIR).glob(pattern))
            if matching_files:
                # Sort by modification time and get the most recent
                src_file = max(matching_files, key=lambda p: p.stat().st_mtime)
                src = str(src_file)
                dst = os.path.join(self.SCRATCHDIR, dst_filename)
                shutil.copy2(src, dst)
                logger.info(f"Copied {src_file.name} -> {dst_filename}")
            else:
                logger.warning(f"No files found matching pattern: {pattern}")
                
        # Handle initial file
        if self.RESTART:
            # Get initial conditions from previous run
            src = os.path.join(self.MSSOUT, self.INIFILE)
            dst = os.path.join(self.SCRATCHDIR, self.INIFILE)
            if os.path.exists(src):
                shutil.copy2(src, dst)
                logger.info(f"Using restart file: {self.INIFILE}")
            else:
                logger.error(f"Restart file not found: {src}")
                return False
        else:
            # Use global data - find the most recent initial file
            ini_pattern = f"{self.MODEL}_ini_mercator_*.nc"
            matching_ini_files = list(Path(self.MSSDIR).glob(ini_pattern))
            if matching_ini_files:
                # Sort by modification time and get the most recent
                src_file = max(matching_ini_files, key=lambda p: p.stat().st_mtime)
                src = str(src_file)
                dst = os.path.join(self.SCRATCHDIR, self.INIFILE)
                shutil.copy2(src, dst)
                logger.info(f"Using initial file from global data: {src_file.name} -> {self.INIFILE}")
            else:
                logger.error(f"No initial files found matching pattern: {ini_pattern}")
                return False
                
        # Copy executable and input files
        files_from_input = [
            self.CODFILE,
            f"{self.MODEL}_forecast.in",
            f"{self.MODEL}_stations.in"
        ]
        
        for filename in files_from_input:
            src = os.path.join(self.INPUTDIR, filename)
            dst = os.path.join(self.SCRATCHDIR, filename)
            if os.path.exists(src):
                shutil.copy2(src, dst)
                if filename == self.CODFILE:
                    os.chmod(dst, 0o755)  # Make executable
                logger.info(f"Copied {filename}")
            else:
                logger.warning(f"File not found: {src}")
                
        return True
        
    def setup_time_management(self):
        """Setup time management in croco_forecast.in"""
        logger.info("Setting up time management...")
        
        # Calculate time parameters
        ndays = self.NDAYS_HIND + self.NDAYS_FCST
        numtimes = ndays * 24 * 3600 // self.DT
        
        numavg = numtimes if self.ND_AVG == -1 else self.ND_AVG * 3600 // self.DT
        numhis = numtimes if self.ND_HIS == -1 else self.ND_HIS * 3600 // self.DT  
        numrst = numtimes if self.ND_RST == -1 else self.ND_RST * 3600 // self.DT
        
        # Read template and substitute values
        template_file = os.path.join(self.INPUTDIR, f"{self.MODEL}_forecast.in")
        output_file = os.path.join(self.SCRATCHDIR, f"{self.MODEL}_forecast.in")
        
        if not os.path.exists(template_file):
            logger.error(f"Template file not found: {template_file}")
            return False
            
        try:
            with open(template_file, 'r') as f:
                content = f.read()
                
            # Substitute placeholders
            replacements = {
                'NUMTIMES': str(numtimes),
                'TIMESTEP': str(self.DT),
                'NFAST': str(self.NFAST),
                'NUMAVG': str(numavg),
                'NUMHIS': str(numhis),
                'NUMRST': str(numrst)
            }
            
            for placeholder, value in replacements.items():
                content = content.replace(placeholder, value)
                
            with open(output_file, 'w') as f:
                f.write(content)
                
            logger.info("Time management setup completed")
            return True
            
        except Exception as e:
            logger.error(f"Error setting up time management: {e}")
            return False
            
    def run_forecast(self):
        """Run the hindcast/forecast simulation"""
        if not self.RUN:
            logger.info("Skipping model run (RUN=0)")
            return True
            
        logger.info("Starting hindcast/forecast run...")
        
        # Check if executable exists
        if not os.path.exists(os.path.join(self.SCRATCHDIR, self.CODFILE)):
            logger.error(f"CROCO executable not found: {self.CODFILE}")
            logger.error("Please ensure the CROCO model is compiled and available")
            return False
        
        # Change to scratch directory
        original_dir = os.getcwd()
        os.chdir(self.SCRATCHDIR)
        
        try:
            # Construct command
            cmd = f"{self.EXEC} {self.CODFILE} {self.MODEL}_forecast.in"
            output_file = f"{self.MODEL}_forecast_{self.datestr}.out"
            
            logger.info(f"Running command: {cmd}")
            logger.info(f"Working directory: {os.getcwd()}")
            logger.info(f"Output will be saved to: {output_file}")
            
            # Run the model
            with open(output_file, 'w') as f:
                # Split command properly for subprocess
                cmd_parts = cmd.split()
                result = subprocess.run(cmd_parts, stdout=f, stderr=subprocess.STDOUT, text=True)
                
            # Check the output file for any immediate errors
            if os.path.exists(output_file):
                with open(output_file, 'r') as f:
                    last_lines = f.readlines()[-10:]  # Read last 10 lines
                    logger.info("Last lines of model output:")
                    for line in last_lines:
                        logger.info(f"  {line.strip()}")
            
            if result.returncode == 0:
                logger.info("Model run completed successfully")
                
                # Store output files
                self.store_output_files()
                return True
            else:
                logger.error(f"Model run failed with return code: {result.returncode}")
                logger.error("Check the output file for detailed error information")
                return False
                
        except Exception as e:
            logger.error(f"Error during model run: {e}")
            return False
        finally:
            os.chdir(original_dir)
            
    def store_output_files(self):
        """Store output forecast files"""
        logger.info("Storing output files...")
        
        # Store restart file for next forecast
        src = os.path.join(self.SCRATCHDIR, self.RSTFILE)
        dst = os.path.join(self.MSSOUT, self.INIFILE)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            logger.info(f"Stored restart file as: {self.INIFILE}")
            
        # Store output files
        output_files = [
            (self.HISFILE, f"{self.MODEL}_his_forecast_{self.datestr}.nc"),
            (self.AVGFILE, f"{self.MODEL}_avg_forecast_{self.datestr}.nc")
        ]
        
        for src_name, dst_name in output_files:
            src = os.path.join(self.SCRATCHDIR, src_name)
            dst = os.path.join(self.MSSOUT, dst_name)
            if os.path.exists(src):
                shutil.copy2(src, dst)
                logger.info(f"Stored: {dst_name}")
            else:
                logger.warning(f"Output file not found: {src}")
                
    def plot_results(self):
        """Generate plots (optional)"""
        if not self.PLOT:
            return True
            
        logger.info("Generating plots...")
        
        # This would need a Python equivalent of plot_forecast_croco.m
        # For now, we'll just log that plotting was requested
        logger.info("Plotting functionality needs Python implementation")
        logger.info("Original MATLAB plotting code needs to be translated")
        
        return True
        
    def generate_tidal_forcing(self):
        """Generate tidal forcing using make_tides.py from CROCO pytools"""
        logger.info("=== Generating tidal forcing ===")
        
        # Check if grid file is available
        if not hasattr(self, 'grdname') or not self.grdname:
            logger.error("Grid file (grdname) not configured")
            return False
            
        if not os.path.exists(self.grdname):
            logger.error(f"Grid file not found: {self.grdname}")
            return False
        
        try:
            # Path to latest CROCO pytools 
            prepro_path = '/DATA/CROCO/croco_pytools-v1.0.3/prepro'
            make_tides_script = os.path.join(prepro_path, 'make_tides.py')
            
            if not os.path.exists(make_tides_script):
                logger.error(f"make_tides.py not found at: {make_tides_script}")
                return False
            
            # Create a custom make_tides configuration file for this forecast
            tides_config_path = os.path.join(self.SCRATCHDIR, 'make_tides_forecast.py')
            
            # Generate configuration content
            config_content = f'''#!/usr/bin/env python3
"""
Auto-generated make_tides.py configuration for CROCO forecast
Generated on: {datetime.datetime.now()}
"""

import netCDF4 as netcdf
import cftime
import numpy as np
import glob
import sys
from os import path
sys.path.append("{prepro_path}/Modules/")
sys.path.append("{prepro_path}/Readers/")
import interp_tools
import croco_class as Croco
import tides_class as Inp

#--- USER CHANGES ---------------------------------------------------------

# Dates
# Initial date
Yini, Mini, Dini = 1980, 1, 1
# Origin year
Yorig, Morig, Dorig = 1980, 1, 1

# Input data information and formating
inputdata = 'tpxo7_croco'
input_dir = '../../DATASETS_CROCOTOOLS/TPXO7/'
input_file = 'TPXO7.nc'
input_type = 'Re_Im'
multi_files = False

# CROCO grid informations
croco_dir = '{os.path.dirname(self.grdname)}/'
croco_grd = '{os.path.basename(self.grdname)}'

# Tide file informations
croco_filename = 'croco_frc.nc'
tides = ['M2','S2','N2','K2','K1','O1','P1','Q1','Mf','Mm']

cur = True  # Set to True if you to compute currents
pot = True  # Set to True if you to compute potiential tides

# Nodal correction
Correction_ssh = True
Correction_uv = True

#--- END USER CHANGES -----------------------------------------------------
'''
            
            # Copy the rest of the original make_tides.py script
            with open(make_tides_script, 'r') as original:
                lines = original.readlines()
                
            # Find where the main script starts (after user changes)
            start_idx = 0
            for i, line in enumerate(lines):
                if '#--- END USER CHANGES' in line:
                    start_idx = i + 1
                    break
            
            # Append the main script logic
            config_content += '\n'.join(lines[start_idx:])
            
            # Write the custom configuration
            with open(tides_config_path, 'w') as f:
                f.write(config_content)
            
            logger.info(f"Created custom tides configuration: {tides_config_path}")
            
            # Change to prepro directory and run the script
            original_dir = os.getcwd()
            original_path = sys.path.copy()
            
            try:
                os.chdir(prepro_path)
                
                # Run the custom tides script
                import subprocess
                result = subprocess.run([
                    sys.executable, tides_config_path
                ], capture_output=True, text=True, cwd=prepro_path)
                
                if result.returncode == 0:
                    logger.info("Tidal forcing generated successfully")
                    logger.info(result.stdout)
                    
                    # Check if the output file was created
                    tides_file = os.path.join(os.path.dirname(self.grdname), 'croco_frc.nc')
                    if os.path.exists(tides_file):
                        logger.info(f"Tidal forcing file created: {tides_file}")
                        return True
                    else:
                        logger.error("Tidal forcing file was not created")
                        return False
                else:
                    logger.error("Error generating tidal forcing:")
                    logger.error(result.stderr)
                    return False
                    
            finally:
                os.chdir(original_dir)
                sys.path = original_path
                
        except Exception as e:
            logger.error(f"Exception in tidal forcing generation: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False
        
    def run_complete_forecast(self):
        """Run the complete forecast workflow"""
        logger.info("=== Starting CROCO Forecast ===")
        logger.info(f"Date: {self.datestr}")
        
        try:
            # Step 1: Clean previous results
            self.clean_previous_results()
            
            # Step 2: Preprocess forcing data
            if not self.preprocess_forcing():
                logger.error("Preprocessing failed")
                return False
                
            # Step 2.5: Generate tidal forcing if enabled
            if self.add_tides_fcst == 1:
                if not self.generate_tidal_forcing():
                    logger.error("Tidal forcing generation failed")
                    return False
                
            # Step 3: Copy files to scratch directory
            if not self.copy_files_to_scratch():
                logger.error("File copying failed")
                return False
                
            # Step 4: Setup time management
            if not self.setup_time_management():
                logger.error("Time management setup failed")
                return False
                
            # Step 5: Run forecast
#            if not self.run_forecast():
#                logger.error("Forecast run failed")
#                return False
                
#            # Step 6: Generate plots
#            if not self.plot_results():
#                logger.error("Plotting failed")
#                return False
                
            logger.info("=== CROCO Forecast Completed Successfully ===")
            return True
            
        except Exception as e:
            logger.error(f"Unexpected error in forecast workflow: {e}")
            return False


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='CROCO Forecast System - Python Version')
    parser.add_argument('--config', help='Configuration file (optional)')
    parser.add_argument('--clean', action='store_true', help='Clean previous results')
    parser.add_argument('--no-preprocess', action='store_true', help='Skip preprocessing')
    parser.add_argument('--restart', action='store_true', help='Restart from previous forecast')
    parser.add_argument('--no-run', action='store_true', help='Skip model run')
    parser.add_argument('--plot', action='store_true', help='Generate plots')
    parser.add_argument('--hindcast-days', type=int, default=1, help='Hindcast days')
    parser.add_argument('--forecast-days', type=int, default=3, help='Forecast days')
    parser.add_argument('--skip-gfs', action='store_true', help='Skip GFS atmospheric forcing and only run OGCM/Mercator processing')
    parser.add_argument('--only-gfs', action='store_true', help='Run only GFS atmospheric forcing and skip OGCM/Mercator processing')

    args = parser.parse_args()
    
    # Create forecast runner
    runner = CROCOForecastRunner(args.config)
    
    # Override defaults with command line arguments
    if args.clean:
        runner.CLEAN = 1
    if args.no_preprocess:
        runner.PRE_PROCESS = 0
    if args.restart:
        runner.RESTART = 1
    if args.no_run:
        runner.RUN = 0
    if args.plot:
        runner.PLOT = 1
    if args.hindcast_days:
        runner.NDAYS_HIND = args.hindcast_days
    if args.forecast_days:
        runner.NDAYS_FCST = args.forecast_days
    # Add skip_gfs and only_gfs flags
    runner.SKIP_GFS = args.skip_gfs
    runner.ONLY_GFS = args.only_gfs
        
    # Run the forecast
    success = runner.run_complete_forecast()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
