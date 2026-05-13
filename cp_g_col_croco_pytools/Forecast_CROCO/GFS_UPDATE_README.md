# GFS Data Access Update

## Summary

The GFS data access has been updated to use the new AWS S3 Open Data method, replacing the deprecated NCEP OpenDAP interface.

## Changes Made

### 1. Updated Data Source
- **Old**: NCEP NOMADS OpenDAP server (`nomads.ncep.noaa.gov`)
- **New**: AWS S3 Open Data bucket (`noaa-gfs-bdp-pds`)

### 2. Updated File Format
- **Old**: NetCDF via OpenDAP (`.nc`)
- **New**: GRIB2 files (`.grib2`)

### 3. Modified Functions

#### `Forecast_tools.py`
- **Updated imports**: Added `boto3`, `botocore`, `xarray`, and `tempfile`
- **Replaced `download_GFS()` function**: Complete rewrite to use AWS S3 API
- **Kept old helper functions**: The OpenDAP-based helper functions (`get_GFS_fname`, `get_GFS_grid`, `get_GFS`, `loaddap`, `readdap`, `getdap`) are still in the file but are no longer used by the main workflow

### 4. New Dependencies Required

Install the following Python packages:

```bash
pip install boto3 xarray cfgrib eccodes
```

Or with conda:

```bash
conda install -c conda-forge boto3 xarray cfgrib eccodes
```

**Note**: `cfgrib` requires the `eccodes` library to read GRIB2 files.

## How It Works

### New Data Access Flow

1. **Connect to AWS S3**: Uses boto3 with unsigned requests (no AWS credentials needed)
2. **Find Latest Forecast**: Searches for the most recent GFS forecast cycle (00, 06, 12, or 18 UTC)
3. **Download GRIB2 Files**: Downloads specific forecast hours as GRIB2 files to a temporary directory
4. **Extract Variables**: Uses xarray with cfgrib engine to read GRIB2 format
5. **Process Data**: Converts to CROCO format (same as before)
6. **Write Output**: Creates the same NetCDF output file as the old method

### Variables Extracted

From GRIB2 files:
- `u10` - 10m u-wind component
- `v10` - 10m v-wind component  
- `t2m` - 2m temperature
- `r2` - 2m relative humidity
- `dswrf` - Downward shortwave radiation flux
- `uswrf` - Upward shortwave radiation flux
- `dlwrf` - Downward longwave radiation flux
- `ulwrf` - Upward longwave radiation flux
- `prate` - Precipitation rate
- `uflx`, `vflx` - Surface momentum fluxes (if available)

### Calculated Variables

- Wind speed: `sqrt(u10^2 + v10^2)`
- Wind stress: Calculated from wind components using bulk formula (Cd = 1.3e-3)
- Net shortwave: `dswrf - uswrf`
- Net longwave: `ulwrf - dlwrf`

## Configuration

No changes to `croco_tools_params.py` are required. The function signature remains the same:

```python
download_GFS(today, lonmin, lonmax, latmin, latmax, FRCST_dir, Yorig, it)
```

## Reference

Based on AWS Open Data samples:
- https://github.com/aws-samples/aws-opendata-samples/blob/main/notebooks/noaa-gfs/noaa_gfs_quickstart.ipynb
- https://registry.opendata.aws/noaa-gfs-bdp-pds/

## Benefits

1. **More Reliable**: AWS S3 is more stable than the deprecated OpenDAP service
2. **No Authentication**: Public bucket requires no credentials
3. **Better Performance**: Direct S3 downloads can be faster
4. **Future-Proof**: AWS Open Data is actively maintained by NOAA

## Backward Compatibility

- Output file format remains the same (NetCDF)
- Variable names and units unchanged
- Integration with existing CROCO forecast system maintained
- The old OpenDAP helper functions are preserved in case they're needed for debugging or reference

## Troubleshooting

### Common Issues

1. **Missing cfgrib**: 
   ```
   pip install cfgrib
   conda install -c conda-forge eccodes
   ```

2. **GRIB2 Read Errors**: Ensure eccodes is properly installed
   ```bash
   conda install -c conda-forge eccodes
   ```

3. **S3 Access Errors**: Check internet connectivity and firewall settings

4. **Variable Not Found**: Some forecast hours may not have all variables. The code includes fallback logic to use previous time step values.

## Testing

To test the updated function:

```python
from datetime import datetime
import Forecast_tools as ft

today = datetime.today()
ft.download_GFS(today, lonmin=-75, lonmax=-60, latmin=35, latmax=45, 
                FRCST_dir='./test_output/', Yorig=2005, it=2)
```

This will download and process GFS data for the specified region and save the output to `./test_output/GFS_XXXXX.nc`.
