#!/usr/bin/env python
import numpy as np
import pyproj


# extent of the voxel space
min_lon = 1.0  # in degrees
max_lon = 1.1
min_lat = 2.3
max_lat = 2.4
min_alt = 0*1000  # in meters
max_alt = 10*1000  # in meters

# resolutions in different directions
resolution_lon = 4
resolution_lat = 4
resolution_alt = 10

input_file_name = "sample_input_data.txt"
input_file_skiprows = 0  # number of header lines to skip in the beginning

# '\s+' is for any number of whitespaces. The r before the quotes is to avoid
# any warnings (in language checkers). This r makes it a raw string.
input_file_delimiter = r"\s+"

# column index in input file. Index start from 0
c_receiver_number = 0  # data is an integer
c_receiver_lon = 1  # data in degrees
c_receiver_lat = 2
c_receiver_alt = 3  # data in meters with respect to reference ellipsoid
c_azimuth = 4  # data in degrees, azimuth and elevation of ray at receiver
c_elevation = 5

# output file names
output_matrix_file_name = "tmp_full_matrix.txt"
validation_file_name = "tmp_validation_matrix.txt"

# plot a particular ray and voxel it has intersected. Index starts from 0.
# If the value is set negative then no ray will be plotted.
plot_epoch_index = 1
# plot_epoch_index = -1

# plot anything at all ? If True, receivers and voxel space top view, 3D view,
# and 3D view with rays will be plotted. Also that ray in plot_epoch_index if
# its value is greater than -1
# master_plot_flag = False
master_plot_flag = True


# corner point values of voxels. These are 1D arrays in increasing order with
# length one greater than respective resolution
lat_points = np.linspace(min_lat, max_lat, resolution_lat + 1)
lon_points = np.linspace(min_lon, max_lon, resolution_lon + 1)
alt_points = np.linspace(min_alt, max_alt, resolution_alt + 1)


# Following block used for calculating radius of earth, from center to
# reference ellipsoid at that particular region (mid point of voxel space).
lon_mid = (min_lon + max_lon)/2
lat_mid = (min_lat + max_lat)/2
# geocentric xyz and lon lat alt frames that will be used for frame transforms
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datumo='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datumo='WGS84')
# getting geocentric xyz of mid point of bottom surface of voxel space
tmpx, tmpy, tmpz = pyproj.transform(lla, ecef,
                                    lon_mid, lat_mid, alt_points[0],
                                    radians=False)
local_earth_rad = np.sqrt(tmpx**2 + tmpy**2 + tmpz**2)  # in meters
# local_earth_rad = 6376972.396554481       # Overwriting if required

x_deg_to_meter_factor = (np.pi/180)*local_earth_rad * np.cos(lat_mid*np.pi/180)
y_deg_to_meter_factor = (np.pi/180)*local_earth_rad
# x_deg_to_meter_factor = 111000  # overwriting
# y_deg_to_meter_factor = 111000


# same as lon, lat_points but in meters starting with 0, i.e. origin at corner.
x_points = (lon_points - min_lon)*x_deg_to_meter_factor
y_points = (lat_points - min_lat)*y_deg_to_meter_factor
# x_points = np.array([0, 5000, 10000, 15000, 20000])  # overwriting
# y_points = np.array([0, 5000, 10000, 15000, 20000])
z_points = alt_points


# This is something that user can decide, but better leave it.
initial_margin = 0.0001  # in meters
initial_delta_margin = initial_margin/10  # adjustments made to margin
delta_margin_factor = 0.1  # multiply factor during delta margin sign change
