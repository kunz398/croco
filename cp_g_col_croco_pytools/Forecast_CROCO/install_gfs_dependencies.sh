#!/bin/bash
# Install required dependencies for updated GFS data access
# Run this script to install boto3, xarray, cfgrib, and eccodes

echo "Installing dependencies for GFS AWS S3 access..."
echo ""

# Check if conda is available
if command -v conda &> /dev/null; then
    echo "Using conda to install dependencies..."
    conda install -c conda-forge boto3 xarray cfgrib eccodes -y
    echo ""
    echo "Dependencies installed via conda!"
else
    echo "Conda not found. Using pip..."
    pip install boto3 xarray cfgrib
    echo ""
    echo "Note: eccodes C library may need to be installed separately."
    echo "On Ubuntu/Debian: sudo apt-get install libeccodes-dev"
    echo "On macOS: brew install eccodes"
    echo "Or use conda: conda install -c conda-forge eccodes"
fi

echo ""
echo "Installation complete!"
echo ""
echo "Testing imports..."
python3 << EOF
try:
    import boto3
    print("✓ boto3 imported successfully")
except ImportError as e:
    print("✗ boto3 import failed:", e)

try:
    import xarray
    print("✓ xarray imported successfully")
except ImportError as e:
    print("✗ xarray import failed:", e)

try:
    import cfgrib
    print("✓ cfgrib imported successfully")
except ImportError as e:
    print("✗ cfgrib import failed:", e)

try:
    from botocore.config import Config
    from botocore import UNSIGNED
    print("✓ botocore imported successfully")
except ImportError as e:
    print("✗ botocore import failed:", e)
EOF

echo ""
echo "If all imports succeeded, you're ready to use the updated GFS download!"
