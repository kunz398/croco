import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Load the netCDF file
file_path = '/DATA/CROCO/croco-v2.1.0/CONFIGS/26-11-2025/SCRATCH/croco_stations.nc'
ds = nc.Dataset(file_path, 'r')

# Read variables
zeta = ds.variables['zeta'][:]  # (ftime, stanum)
scrum_time = ds.variables['scrum_time'][:]  # time in seconds
lat = ds.variables['lat'][:]
lon = ds.variables['lon'][:]

# Convert time to days (or use datetime if you have a reference date)
time_days = scrum_time / 86400.0  # Convert seconds to days

# Create subplots (3x3 grid for 9 stations)
fig, axes = plt.subplots(3, 3, figsize=(15, 12))
axes = axes.flatten()

# Plot each station
for i in range(9):
    ax = axes[i]
    
    # Plot zeta time series
    ax.plot(time_days, zeta[:, i], 'b-', linewidth=0.8)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Labels and title
    ax.set_title(f'Station {i+1}', fontsize=11, fontweight='bold')
    ax.set_xlabel('Time (days)', fontsize=9)
    ax.set_ylabel('Zeta (m)', fontsize=9)
    
    # Add lat/lon text in top-right corner
    textstr = f'Lat: {lat[i]:.3f}°\nLon: {lon[i]:.3f}°'
    ax.text(0.98, 0.97, textstr, transform=ax.transAxes,
            fontsize=8, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    # Adjust tick label size
    ax.tick_params(labelsize=8)

# Adjust layout
plt.suptitle('Sea Level (Zeta) Time Series at CROCO Stations', 
             fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout()

# Save figure
plt.savefig('croco_zeta_timeseries.png', dpi=300, bbox_inches='tight')
print("Figure saved as 'croco_zeta_timeseries.png'")

# Show plot
plt.show()

# Close dataset
ds.close()