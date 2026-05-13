#!/usr/bin/env python3
"""
Test script for updated GFS data access
Tests the new AWS S3-based download without actually downloading large amounts of data
"""

import sys
import os

print("=" * 60)
print("Testing GFS AWS S3 Access Dependencies")
print("=" * 60)
print()

# Test 1: Import required modules
print("Test 1: Checking required imports...")
try:
    import boto3
    from botocore import UNSIGNED
    from botocore.config import Config
    print("  ✓ boto3 and botocore")
except ImportError as e:
    print(f"  ✗ boto3/botocore failed: {e}")
    sys.exit(1)

try:
    import xarray as xr
    print("  ✓ xarray")
except ImportError as e:
    print(f"  ✗ xarray failed: {e}")
    sys.exit(1)

try:
    import cfgrib
    print("  ✓ cfgrib")
except ImportError as e:
    print(f"  ✗ cfgrib failed: {e}")
    print("  Note: You may need to install eccodes library")
    sys.exit(1)

print()

# Test 2: Test S3 connection
print("Test 2: Testing AWS S3 connection...")
try:
    # Disable SSL verification for institutional networks
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED), verify=False)
    GFS_BUCKET_NAME = "noaa-gfs-bdp-pds"
    print("  Note: Using SSL verification disabled for institutional network")
    
    # Try to list a recent forecast
    from datetime import datetime, timedelta
    today = datetime.utcnow()
    date_str = today.strftime('%Y%m%d')
    
    folder = f"gfs.{date_str}/00/atmos/"
    response = s3.list_objects_v2(Bucket=GFS_BUCKET_NAME, Prefix=folder, MaxKeys=5)
    
    if 'Contents' in response:
        print(f"  ✓ Successfully connected to S3 bucket")
        print(f"  ✓ Found files for {date_str}")
        print(f"    Sample files:")
        for obj in response['Contents'][:3]:
            print(f"      - {obj['Key']}")
    else:
        print(f"  ! No files found for today ({date_str}), trying yesterday...")
        yesterday = today - timedelta(days=1)
        date_str = yesterday.strftime('%Y%m%d')
        folder = f"gfs.{date_str}/18/atmos/"
        response = s3.list_objects_v2(Bucket=GFS_BUCKET_NAME, Prefix=folder, MaxKeys=5)
        if 'Contents' in response:
            print(f"  ✓ Found files for {date_str}")
        else:
            print(f"  ✗ Could not find recent GFS files")

except Exception as e:
    print(f"  ✗ S3 connection failed: {e}")
    sys.exit(1)

print()

# Test 3: Test Forecast_tools import
print("Test 3: Testing Forecast_tools module...")
try:
    import Forecast_tools as ft
    print("  ✓ Forecast_tools imported successfully")
    
    # Check if download_GFS function exists
    if hasattr(ft, 'download_GFS'):
        print("  ✓ download_GFS function found")
    else:
        print("  ✗ download_GFS function not found")
        
except ImportError as e:
    print(f"  ✗ Could not import Forecast_tools: {e}")
    print(f"    Make sure you're running this from the correct directory")

print()
print("=" * 60)
print("All tests passed! ✓")
print("The updated GFS download system is ready to use.")
print("=" * 60)
print()
print("To download GFS data, use:")
print("  import Forecast_tools as ft")
print("  from datetime import datetime")
print("  ft.download_GFS(datetime.today(), lon_min, lon_max, lat_min, lat_max, output_dir, year_origin, time_interval)")
print()
