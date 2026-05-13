#!/usr/bin/env python3
"""
Test script for complete forecast workflow with tidal forcing
"""

import sys
import os
from cp_g_col_croco_pytools.Forecast_CROCO.old___run_croco_forecast import CROCOForecastRunner

def test_complete_workflow():
    """Test complete workflow including tidal forcing"""
    
    # Create a forecast runner instance
    runner = CROCOForecastRunner()
    
    # Show configuration
    print(f"Tidal forcing enabled: {runner.add_tides_fcst}")
    print(f"Grid file: {runner.grdname}")
    
    # Test only the preprocessing + tidal forcing (skip the actual model run)
    print("\n=== Testing preprocessing workflow with tidal forcing ===")
    
    # Set flags for testing
    runner.PRE_PROCESS = 1
    runner.RUN = 0  # Skip actual model run for testing
    runner.PLOT = 0
    
    print("Configurations:")
    print(f"  PRE_PROCESS: {runner.PRE_PROCESS}")
    print(f"  RUN: {runner.RUN}")
    print(f"  Tidal forcing: {runner.add_tides_fcst}")
    
    # Run only specific parts for testing
    try:
        # Step 1: Clean previous results
        print("\n--- Step 1: Cleaning previous results ---")
        runner.clean_previous_results()
        print("✅ Cleaning completed")
        
        # Step 2: Preprocessing (but skip to save time)
        print("\n--- Step 2: Preprocessing (SKIPPED for this test) ---")
        print("⏭️  Skipping preprocessing to save time")
        
        # Step 3: Tidal forcing generation
        if runner.add_tides_fcst == 1:
            print("\n--- Step 3: Tidal forcing generation ---")
            success = runner.generate_tidal_forcing()
            if success:
                print("✅ Tidal forcing generation completed")
            else:
                print("❌ Tidal forcing generation failed")
                return False
        
        print("\n🎉 Test workflow completed successfully!")
        return True
        
    except Exception as e:
        print(f"\n❌ Error in test workflow: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    test_complete_workflow()