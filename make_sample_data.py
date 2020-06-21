import numpy as np

number_of_receivers = 20
number_of_rays = 100

receiver_min_lon = 1.0  # in degrees
receiver_max_lon = 1.1
receiver_min_lat = 2.3
receiver_max_lat = 2.4
receiver_min_alt = -1000  # in meters
receiver_max_alt = 2000

save_file_name = 'sample_input_data.txt'

min_azimuth = 0  # deg
max_azimuth = 360
min_elevation = 80
max_elevation = 90

receiver_number_list = np.arange(number_of_receivers) + 1
receiver_lon_list = np.random.uniform(receiver_min_lon, receiver_max_lon,
                                      number_of_receivers)
receiver_lat_list = np.random.uniform(receiver_min_lat, receiver_max_lat,
                                      number_of_receivers)
receiver_alt_list = np.random.uniform(receiver_min_alt, receiver_max_alt,
                                      number_of_receivers)

data_order_index = np.random.choice(np.arange(number_of_receivers),
                                    number_of_rays)
receiver_number = receiver_number_list[data_order_index]
receiver_lon = receiver_lon_list[data_order_index]
receiver_lat = receiver_lat_list[data_order_index]
receiver_alt = receiver_alt_list[data_order_index]
azimuth = np.random.uniform(min_azimuth, max_azimuth, number_of_rays)
elevation = np.arcsin(np.random.uniform(
            np.sin(min_elevation*np.pi/180), np.sin(max_elevation*np.pi/180),
            number_of_rays)) * 180/np.pi

data = np.column_stack((receiver_number, receiver_lon))
data = np.column_stack((data, receiver_lat))

data = np.column_stack((
    receiver_number,
    receiver_lon,
    receiver_lat,
    receiver_alt,
    azimuth,
    elevation
    ))

np.savetxt(save_file_name, data, delimiter=' ')
