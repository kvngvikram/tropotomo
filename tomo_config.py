#!/usr/bin/env python
import numpy as np
import pyproj

ray_cutoff_height = 8000   # geodetic height in m

min_lon = 1.0  # in degrees
max_lon = 1.1
min_lat = 2.3
max_lat = 2.4
min_alt = 0*1000  # in meters
max_alt = 10*1000  # in meters
resolution_lon = 4
resolution_lat = 4
resolution_alt = 10

lat_points = np.linspace(min_lat, max_lat, resolution_lat + 1)
lon_points = np.linspace(min_lon, max_lon, resolution_lon + 1)
alt_points = np.linspace(min_alt, max_alt, resolution_alt + 1)

# | geocentric xyz and lon lat alt frames
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datumo='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datumo='WGS84')


# | function that converts lon lat alt to geocentric xyz
def get_ecef(lon, lat, height):
    return pyproj.transform(lla, ecef, lon, lat, height, radians=False)


lon_mid = (min_lon + max_lon)/2
lat_mid = (min_lat + max_lat)/2
tmpx, tmpy, tmpz = get_ecef(
        lon_mid,
        lat_mid,
        alt_points[0])
earth_rad = np.sqrt(tmpx**2 + tmpy**2 + tmpz**2)
# earth_rad = 6376972.396554481       # Overwriting R

x_deg_to_meter_factor = (np.pi/180)*earth_rad * np.cos(lat_mid * np.pi/180)
y_deg_to_meter_factor = (np.pi/180)*earth_rad
# x_deg_to_meter_factor = 111000
# y_deg_to_meter_factor = 111000


x_points = (lon_points - min_lon)*x_deg_to_meter_factor
y_points = (lat_points - min_lat)*y_deg_to_meter_factor
# x_points = np.array([0, 5000, 10000, 15000, 20000])
# y_points = np.array([0, 5000, 10000, 15000, 20000])
z_points = alt_points


input_file_name = "sample_input_data.txt"
input_file_skiprows = 1
# '\s+' is for any number of whitespaces. The r before the quotes is to avoid
# any warnings (in language checkers). This r makes it a raw string.
input_file_delimiter = r"\s+"
c_receiver_number = 0
c_receiver_lon = 1
c_receiver_lat = 2
c_receiver_alt = 3
c_azimuth = 4
c_elevation = 5

output_matrix_file_name = "tmp_full_matrix.txt"
validation_file_name = "tmp_validation_matrix.txt"

# plot this particular ray and voxel it has intersected. Index starts from 0.
# If the value is set negative then no ray will be plotted.
plot_epoch_index = 10
# plot_epoch_index = -1

# plot anything at all ?
# master_plot_flag = False
master_plot_flag = True

# This is something that user can decide, but better leave it.
initial_margin = 0.0001  # in meters
initial_delta_margin = initial_margin/10  # adjustments made to margin
delta_margin_factor = 0.1  # multiply factor during delta margin sign change
