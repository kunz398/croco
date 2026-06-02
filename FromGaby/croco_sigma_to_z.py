#!/usr/bin/env python3
"""
Python script to convert CROCO sigma-coordinate output to z-levels using xcroco.
Reproduces the functionality of format_croco_z.sh
"""

import sys
import os
import logging
from pathlib import Path

def _ensure_xcroco_on_path():
    """Ensure xcroco is importable from install, env override, or known local roots."""
    try:
        import xcroco  # noqa: F401
        return
    except ModuleNotFoundError:
        pass

    script_dir = Path(__file__).resolve().parent
    env_path = os.environ.get('XCROCO_PATH', '').strip()

    candidates = [
        Path(env_path) if env_path else None,
        script_dir.parent / 'cp_g_col_croco_pytools',
        Path('/DATA/CROCO/cp_g_col_croco_pytools'),
    ]

    for candidate in candidates:
        if candidate is None:
            continue
        pkg_root = candidate / 'xcroco' / '__init__.py'
        if pkg_root.exists():
            candidate_str = str(candidate)
            if candidate_str not in sys.path:
                sys.path.insert(0, candidate_str)
            return


_ensure_xcroco_on_path()

import xarray as xr
import numpy as np

# Logger — writes to stdout so it appears in run_forecast.log via the container pipe
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Import xcroco modules
try:
    from xcroco.model import Model
    import xcroco.inout as io
    import xcroco.gridop as gop
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "xcroco is not available. Install it or set XCROCO_PATH to a folder containing xcroco/."
    ) from exc


def _get_rho_mask(ds, grid_file):
    """Return CROCO rho-point land/sea mask as boolean ocean mask (True=ocean)."""
    candidate_names = ['mask_rho', 'mask', 'rmask']

    def _find_mask(obj):
        for name in candidate_names:
            if name in obj:
                return obj[name]
        return None

    mask = _find_mask(ds.data_vars)
    if mask is None:
        mask = _find_mask(ds.coords)

    if mask is None:
        with xr.open_dataset(grid_file) as grd:
            mask = _find_mask(grd.variables)
            if mask is not None:
                mask = mask.load()

    if mask is None:
        return None

    rename_map = {}
    if 'eta_rho' in mask.dims:
        rename_map['eta_rho'] = 'y'
    if 'xi_rho' in mask.dims:
        rename_map['xi_rho'] = 'x'
    if rename_map:
        mask = mask.rename(rename_map)

    keep_dims = [dim for dim in mask.dims if dim in ['y', 'x']]
    if len(keep_dims) < 2:
        return None

    drop_dims = [dim for dim in mask.dims if dim not in ['y', 'x']]
    if drop_dims:
        mask = mask.isel({dim: 0 for dim in drop_dims})

    return mask > 0.5


def croco_sigma_to_z(input_file, output_file, grid_file, 
                     z_levels=None, chunks={'t': 1}, variables=None):
    """
    Convert CROCO sigma-coordinate output to z-levels.
    
    Parameters
    ----------
    input_file : str
        Path to input CROCO history file (croco_his.nc)
    output_file : str
        Path to output NetCDF file with z-level interpolation
    grid_file : str
        Path to CROCO grid file (croco_grd.nc)
    z_levels : array-like, optional
        Array of z-levels in meters (negative for depth below surface).
        Default: [-5, -10, -20, -30, -50, -100, -300, -500, -1000]
    chunks : dict, optional
        Chunk specification for dask. Default: {'t': 1}
    variables : list, optional
        List of variables to interpolate. Default: ['temp', 'salt', 'xcur', 'ycur']
        
    Returns
    -------
    ds_out : xarray.Dataset
        Output dataset with interpolated variables on z-levels
        
    Examples
    --------
    Basic usage with default z-levels:
    >>> ds = croco_sigma_to_z('croco_his.nc', 'output_z.nc', 'croco_grd.nc')
    
    Custom z-levels:
    >>> z_levels = [-10, -20, -50, -100, -200]
    >>> ds = croco_sigma_to_z('croco_his.nc', 'output_z.nc', 'croco_grd.nc', 
    ...                       z_levels=z_levels)
    
    Process large files with custom chunking:
    >>> ds = croco_sigma_to_z('croco_his.nc', 'output_z.nc', 'croco_grd.nc',
    ...                       chunks={'t': 1})
    """
    
    # Default z-levels if not provided
    if z_levels is None:
        z_levels = np.array([-5, -10, -20, -30, -50, -100, -300, -500, -1000])
    else:
        z_levels = np.array(z_levels)
    
    # Default variables if not provided
    if variables is None:
        variables = ['temp', 'salt', 'xcur', 'ycur']
    
    logger.info("=== sigma-to-z post-processing started ===")
    logger.info("Setting up CROCO model...")
    # Create model instance (use croco_native for standard CROCO output)
    croco = Model("croco_native")
    
    logger.info(f"Loading CROCO output file: {input_file}")
    # Open files using xcroco's io.open_files
    # Note: io.open_files returns (ds, grid) tuple
    # xcur and ycur are the u and v velocity components in xcroco
    [ds, xgrid] = io.open_files(croco, grid_file, [input_file], 
                           chunks=chunks,  # Use 't' not 'time' for CROCO native files
                           drop_variables=[])  

    logger.info(f"Variables in dataset: {list(ds.data_vars.keys())}")
    logger.info(f"Coordinates: {list(ds.coords.keys())}")

    # Get z coordinates at rho points
    z_rho = gop.get_z(croco, ds=ds, z_sfc=ds.z_sfc, xgrid=xgrid, hgrid='r', vgrid='r')
    logger.info(f"z_rho dims: {z_rho.dims}")
    logger.info(f"z_rho min: {z_rho.min().values:.2f}, max: {z_rho.max().values:.2f}, NaN%: {100*np.isnan(z_rho.values).sum()/z_rho.values.size:.1f}%")
    logger.info(f"z_rho upper layer min: {z_rho[:,49,:,:].min().values:.2f}, max: {z_rho[:,49,:,:].max().values:.2f}")

    logger.info("Interpolating temperature to z-levels...")
    temp_z = gop.isoslice(ds.temp, z_levels, xgrid, target_data=z_rho, axis="z")
    logger.info(f"temp_z: {temp_z.dims}, {temp_z.shape}")

    logger.info("Interpolating salinity to z-levels...")
    salt_z = gop.isoslice(ds.salt, z_levels, xgrid, target_data=z_rho, axis="z")
    logger.info(f"salt_z: {salt_z.dims}, {salt_z.shape}")

    # Get z coordinates at u points for u velocity
    logger.info("Getting z coordinates at u points...")
    z_u = gop.get_z(croco, ds=ds, z_sfc=ds.z_sfc, xgrid=xgrid, hgrid='u', vgrid='r')
    logger.info(f"z_u dims: {z_u.dims}")

    logger.info("Interpolating u velocity to z-levels...")
    u_z = gop.isoslice(ds.xcur, z_levels, xgrid, target_data=z_u, axis="z")
    logger.info(f"u_z: {u_z.dims}, {u_z.shape}")

    # Get z coordinates at v points for v velocity
    logger.info("Getting z coordinates at v points...")
    z_v = gop.get_z(croco, ds=ds, z_sfc=ds.z_sfc, xgrid=xgrid, hgrid='v', vgrid='r')
    logger.info(f"z_v dims: {z_v.dims}")

    logger.info("Interpolating v velocity to z-levels...")
    v_z = gop.isoslice(ds.ycur, z_levels, xgrid, target_data=z_v, axis="z")
    logger.info(f"v_z: {v_z.dims}, {v_z.shape}")

    # Now interpolate u and v from their staggered grids to rho points
    logger.info("Interpolating u from u-grid to rho-grid...")
    u_z_rho = xgrid.interp(u_z, 'x', boundary='extend')
    # Rename z_u to z_r since it's now on rho grid
    u_z_rho = u_z_rho.rename({'z_u': 'z_r'})
    # Copy attributes from original xcur
    u_z_rho.attrs = ds.xcur.attrs.copy()
    logger.info(f"u_z_rho: {u_z_rho.dims}, {u_z_rho.shape}")

    logger.info("Interpolating v from v-grid to rho-grid...")
    v_z_rho = xgrid.interp(v_z, 'y', boundary='extend')
    # Rename z_v to z_r since it's now on rho grid
    v_z_rho = v_z_rho.rename({'z_v': 'z_r'})
    # Copy attributes from original ycur
    v_z_rho.attrs = ds.ycur.attrs.copy()
    logger.info(f"v_z_rho: {v_z_rho.dims}, {v_z_rho.shape}")

    # Apply rho-point land mask after interpolation to prevent values over land
    rho_ocean_mask = _get_rho_mask(ds, grid_file)
    if rho_ocean_mask is not None:
        u_z_rho = u_z_rho.where(rho_ocean_mask)
        v_z_rho = v_z_rho.where(rho_ocean_mask)
        logger.info("Applied rho-point land mask to u and v.")
    else:
        logger.warning("rho-point mask not found; u/v land masking was not applied.")

    # Create a new dataset with interpolated variables
    ds_out = xr.Dataset()
    ds_out['temperature'] = temp_z
    ds_out['salinity'] = salt_z
    ds_out['u'] = u_z_rho
    ds_out['v'] = v_z_rho

    # Rename dimensions to CF-compliant names
    ds_out = ds_out.rename({
        't': 'time',
        'z_r': 'depth'
    })

    # The 'depth' dimension now exists; assign values and attributes to the depth coordinate
    ds_out['depth'] = z_levels.astype('float64')
    ds_out['depth'].attrs = {
        'long_name': 'depth below sea surface',
        'standard_name': 'depth',
        'units': 'm',
        'positive': 'up',
        'axis': 'Z',
        'comment': 'negative values indicate depth below sea surface'
    }

    # Ensure x and y coordinates are float64 (convert dimension indices to floats)
    if 'x' in ds_out.coords:
        ds_out['x'] = ds_out['x'].astype('float64')
    if 'y' in ds_out.coords:
        ds_out['y'] = ds_out['y'].astype('float64')

    # Copy time coordinate attributes from original dataset
    if 't' in ds.coords:
        # Assign the time coordinate values and attributes
        ds_out['time'] = ds.t.values
        ds_out['time'].attrs = ds['t'].attrs.copy()
        # Remove _FillValue if it exists
        if '_FillValue' in ds_out['time'].attrs:
            del ds_out['time'].attrs['_FillValue']
        ds_out['time'].attrs.update({
            'units': 'seconds since 2005-01-01 00:00:00',
            'calendar': 'gregorian',
            'axis': 'T',
            'standard_name': 'time'
        })

    # Copy spatial coordinates (lat/lon) from original dataset and rename to latitude/longitude
    if 'lat' in ds.coords:
        lat_clean = ds.coords['lat']
        ds_out.coords['latitude'] = lat_clean
        ds_out['latitude'].attrs = ds['lat'].attrs.copy()
        ds_out['latitude'].attrs['standard_name'] = 'latitude'
        # Remove _FillValue if it exists
        if '_FillValue' in ds_out['latitude'].attrs:
            del ds_out['latitude'].attrs['_FillValue']

    if 'lon' in ds.coords:
        lon_clean = ds.coords['lon']
        ds_out.coords['longitude'] = lon_clean
        ds_out['longitude'].attrs = ds['lon'].attrs.copy()
        ds_out['longitude'].attrs['standard_name'] = 'longitude'
        # Remove _FillValue if it exists
        if '_FillValue' in ds_out['longitude'].attrs:
            del ds_out['longitude'].attrs['_FillValue']

    # Drop the old 'lat' and 'lon' coordinates now that we have 'latitude' and 'longitude'
    if 'lat' in ds_out.coords:
        ds_out = ds_out.drop_vars('lat')
    if 'lon' in ds_out.coords:
        ds_out = ds_out.drop_vars('lon')
     
    # Update variable attributes for CF compliance
    # Temperature
    ds_out['temperature'].attrs.update({
        'standard_name': 'sea_water_potential_temperature',
        'long_name': 'potential temperature',
        'coordinates': 'time depth latitude longitude',
        'units': 'degree_Celsius',
        '_FillValue': np.nan
    })

    # Salinity
    ds_out['salinity'].attrs.update({
        'standard_name': 'sea_water_salinity',
        'long_name': 'salinity',
        'coordinates': 'time depth latitude longitude',
        'units': '0.001', #CF convention uses dimensionless "1" for salinity (or "0.001" for practical salinity)
        '_FillValue': np.nan
    })

    # U velocity
    ds_out['u'].attrs.update({
        'standard_name': 'eastward_sea_water_velocity',
        'long_name': 'u-momentum component at rho points',
        'coordinates': 'time depth latitude longitude',
        'units': 'm s-1',
        'comment': 'interpolated from u-grid to rho-grid',
        '_FillValue': np.nan
    })

    # V velocity
    ds_out['v'].attrs.update({
        'standard_name': 'northward_sea_water_velocity',
        'long_name': 'v-momentum component at rho points',
        'coordinates': 'time depth latitude longitude',
        'units': 'm s-1',
        'comment': 'interpolated from v-grid to rho-grid',
        '_FillValue': np.nan
    })

    # Global attributes for CF compliance
    ds_out.attrs = {
        'Conventions': 'CF-1.8',
        'title': 'CROCO model output interpolated to z-levels',
        'institution': ds.attrs.get('institution', 'Unknown'),
        'source': 'CROCO (Coastal and Regional Ocean COmmunity model) configured and ran by Gaby Mayorga-Adame for Pacific Community (SPC)',
        'history': f'Interpolated from sigma to z-levels on {np.datetime64("today")}',
        'references': 'https://www.croco-ocean.org',
        'comment': f'Interpolated to {len(z_levels)} z-levels from original CROCO sigma coordinates',
        'vertical_interpolation_method': 'Linear interpolation using xcroco gop.isoslice',
        'horizontal_interpolation_method': 'Linear interpolation using xgcm grid.interp for u/v to rho points'
    }

    # Copy some important attributes from original dataset
    for attr in ['type', 'title']:
        if attr in ds.attrs:
            ds_out.attrs[f'original_{attr}'] = ds.attrs[attr]

    # Drop the old 'z' coordinate from the dataset if it still exists
    if 'z' in ds_out.coords:
        ds_out = ds_out.drop_vars('z')

    # Add zeta (free-surface elevation) — already on rho grid, no interpolation needed.
    # xcroco may not expose zeta in ds.data_vars, so read it straight from the raw file.
    zeta_added = False
    try:
        with xr.open_dataset(input_file) as ds_raw:
            if 'zeta' in ds_raw.data_vars:
                zeta = ds_raw['zeta'].load()
                # Rename CROCO native dims to match ds_out: eta_rho→y, xi_rho→x, ocean_time→time
                rename_map = {}
                if 'eta_rho' in zeta.dims:
                    rename_map['eta_rho'] = 'y'
                if 'xi_rho' in zeta.dims:
                    rename_map['xi_rho'] = 'x'
                # CROCO time dimension can be 'ocean_time' or 'time' in the raw file
                for tdim in ('ocean_time', 'scrum_time', 'time'):
                    if tdim in zeta.dims:
                        rename_map[tdim] = 'time'
                        break
                if rename_map:
                    zeta = zeta.rename(rename_map)
                # Replace time coordinate values with the CF-compliant ones already in ds_out
                zeta['time'] = ds_out['time']
                # Attach spatial auxiliary coordinates so THREDDS treats it as geo2d
                if 'latitude' in ds_out.coords:
                    zeta = zeta.assign_coords(latitude=ds_out['latitude'])
                if 'longitude' in ds_out.coords:
                    zeta = zeta.assign_coords(longitude=ds_out['longitude'])
                zeta.attrs.update({
                    'standard_name': 'sea_surface_height_above_geoid',
                    'long_name': 'free-surface elevation',
                    'units': 'm',
                })
                # Do NOT set 'coordinates' in attrs — xarray encodes auxiliary
                # coordinates automatically and raises ValueError if it's also in attrs.
                ds_out['zeta'] = zeta
                zeta_added = True
                logger.info(f"Added zeta to output dataset with dims {zeta.dims}")
            else:
                logger.warning("zeta not found in raw history file; skipping.")
    except Exception as e:
        logger.warning(f"Could not read zeta from raw file: {e}")

    # Set encoding to ensure coordinates are proper types and have no _FillValue
    encoding = {
        'x': {'dtype': 'float64', '_FillValue': None},
        'y': {'dtype': 'float64', '_FillValue': None},
        'depth': {'dtype': 'float64', '_FillValue': None},
        'time': {'_FillValue': None},
        'latitude': {'_FillValue': None},
        'longitude': {'_FillValue': None},
    }
    if zeta_added:
        encoding['zeta'] = {'_FillValue': float('nan')}

    # Write to output file
    if os.path.exists(output_file):
        os.remove(output_file)
    
    ds_out.to_netcdf(output_file, compute=True, encoding=encoding)
    logger.info(f"Saved {output_file}")
    logger.info("=== sigma-to-z post-processing complete ===")
    
    return ds_out


if __name__ == '__main__':
    import datetime
    script_dir = Path(__file__).parent.resolve()
    today = datetime.date.today().strftime('%d-%m-%Y')
    scratch_dir = script_dir.parent / today / 'SCRATCH'

    input_file = str(scratch_dir / 'croco_his.nc')
    grid_file = str(scratch_dir / 'croco_grd.nc')
    output_file = str(scratch_dir / 'temp_salt_uv_z.nc')

    # Original hardcoded paths (kept for reference)
    # gpath = '/media/claudiaa/Elements/eDRIVE/02-02-2026/SCRATCH/'
    # path = '/DATA/CROCO/FromDMZ/'
    # input_file = os.path.join(path, 'croco_his.00000.nc')
    # grid_file = os.path.join(gpath, 'croco_grd.nc')
    # output_file = os.path.join(path, '1days_temp_salt_uv_z.nc')

    logger.info(f"Input:  {input_file}")
    logger.info(f"Grid:   {grid_file}")
    logger.info(f"Output: {output_file}")

    # Run the function
    ds_out = croco_sigma_to_z(input_file, output_file, grid_file)
    logger.info(f"croco_sigma_to_z completed")
