#!/usr/bin/env python3
"""
Test script to manually generate tidal forcing and check file creation
"""

import sys
import os
from cp_g_col_croco_pytools.Forecast_CROCO.old___run_croco_forecast import CROCOForecastRunner

def test_tidal_manual():
    """Test manual tidal forcing generation"""
    
    # Create a forecast runner instance
    runner = CROCOForecastRunner()
    
    print(f"Configuration:")
    print(f"  add_tides_fcst: {runner.add_tides_fcst}")
    print(f"  Grid file: {runner.grdname}")
    print(f"  CROCO_FILES dir: {runner.MSSDIR}")
    
    if runner.add_tides_fcst == 1:
        print("\n=== Manual Tidal Forcing Generation ===")
        
        # Test the tidal forcing generation
        success = runner.generate_tidal_forcing()
        
        if success:
            print("✅ Tidal forcing generation completed!")
            
            # Check if file was created
            tidal_file = os.path.join(runner.MSSDIR, 'croco_frc.nc')
            if os.path.exists(tidal_file):
                print(f"✅ Tidal file found: {tidal_file}")
                import subprocess
                result = subprocess.run(['ls', '-lh', tidal_file], capture_output=True, text=True)
                print(f"File size: {result.stdout.strip()}")
            else:
                print(f"❌ Tidal file NOT found at: {tidal_file}")
                
                # Check if it was created elsewhere
                possible_locations = [
                    os.path.dirname(runner.grdname),
                    runner.SCRATCHDIR,
                    os.getcwd()
                ]
                
                for location in possible_locations:
                    test_path = os.path.join(location, 'croco_frc.nc')
                    if os.path.exists(test_path):
                        print(f"🔍 Found tidal file at: {test_path}")
                        import subprocess
                        result = subprocess.run(['ls', '-lh', test_path], capture_output=True, text=True)
                        print(f"File size: {result.stdout.strip()}")
                        break
                else:
                    print("❌ Tidal file not found in any expected location")
        else:
            print("❌ Tidal forcing generation failed!")
    else:
        print("Tidal forcing is disabled")

if __name__ == "__main__":
    test_tidal_manual()