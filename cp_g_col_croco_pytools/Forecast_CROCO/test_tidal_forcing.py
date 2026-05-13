#!/usr/bin/env python3
"""
Test script for tidal forcing generation
"""

import sys
import os
from cp_g_col_croco_pytools.Forecast_CROCO.old___run_croco_forecast import CROCOForecastRunner

def test_tidal_forcing():
    """Test tidal forcing generation"""
    
    # Create a forecast runner instance
    runner = CROCOForecastRunner()
    
    # Show configuration
    print(f"Tidal forcing enabled: {runner.add_tides_fcst}")
    print(f"Grid file: {runner.grdname}")
    print(f"Grid file exists: {os.path.exists(runner.grdname)}")
    
    if runner.add_tides_fcst == 1:
        print("\n=== Testing tidal forcing generation ===")
        
        # Test the tidal forcing generation
        success = runner.generate_tidal_forcing()
        
        if success:
            print("✅ Tidal forcing generation completed successfully!")
        else:
            print("❌ Tidal forcing generation failed!")
            return False
    else:
        print("Tidal forcing is disabled in configuration")
        
    return True

if __name__ == "__main__":
    test_tidal_forcing()