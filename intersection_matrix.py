#!/usr/bin/env python
import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.style import use as mpluse

# get variables from config file
from tomo_config import (
        # file names
        input_file_name, input_file_skiprows, input_file_delimiter,
        output_matrix_file_name, validation_file_name,

        # column numbers in input data file
        c_receiver_lon, c_receiver_lat, c_receiver_alt, c_azimuth, c_elevation,
        c_receiver_number,

        # conversion factors from lat lon to meters in that area
        x_deg_to_meter_factor, y_deg_to_meter_factor,

        # voxel space parameters
        min_lon, min_lat, x_points, y_points, z_points,

        ray_cutoff_height,

        # index of ray that should be plotted and should stuff be plotted?
        plot_epoch_index, master_plot_flag,

        # margin parameters used during intersection length calculations
        initial_margin, initial_delta_margin, delta_margin_factor
        )

plot_epoch_index = -1 if master_plot_flag else plot_epoch_index

# --------------------------------------------------------
# | data reading and assignement
# --------------------------------------------------------
data = np.asarray(read_csv(
        input_file_name, delimiter=input_file_delimiter, header=None,
        skiprows=input_file_skiprows
        ))

receiver_number = data[:, c_receiver_number]
receiver_lon = data[:, c_receiver_lon]
receiver_lat = data[:, c_receiver_lat]
receiver_alt = data[:, c_receiver_alt]

azimuth = data[:, c_azimuth]
elevation = data[:, c_elevation]

del data  # this data is not required from now


# --------------------------------------------------------
# | line equation calculations
# --------------------------------------------------------
x0 = (receiver_lon - min_lon)*x_deg_to_meter_factor
y0 = (receiver_lat - min_lat)*y_deg_to_meter_factor
z0 = receiver_alt
sr = 10
x1 = sr*np.cos(elevation * np.pi/180)*np.sin(azimuth * np.pi/180)
y1 = sr*np.cos(elevation * np.pi/180)*np.cos(azimuth * np.pi/180)
z1 = sr*np.sin(elevation * np.pi/180)


# | if plotting is required for a particular ray and its voxel intersections
if(plot_epoch_index >= 0):
    fig = plt.figure('ray intersections')
    ax = plt.axes(projection=Axes3D.name)
    tmpx = x_points
    tmpy = y_points
    tmpX, tmpY = np.meshgrid(tmpx, tmpy)
    tmpZ = tmpX*0 + z_points.min()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')


# --------------------------------------------------------
# | calculating intersection points for all combinations
# --------------------------------------------------------
x_alpha = (x_points[None, :] - x0[:, None])/x1[:, None]
y_alpha = (y_points[None, :] - y0[:, None])/y1[:, None]
z_alpha = (z_points[None, :] - z0[:, None])/z1[:, None]
x_alpha[np.isinf(x_alpha)] = -1
y_alpha[np.isinf(y_alpha)] = -1
z_alpha[np.isinf(z_alpha)] = -1
x_alpha[np.isnan(x_alpha)] = -1
y_alpha[np.isnan(y_alpha)] = -1
z_alpha[np.isnan(z_alpha)] = -1

# --------------------------------------------------------
# | intersection length calculations
# --------------------------------------------------------
x_index = np.linspace(0, x_points.size-2, x_points.size-1, dtype=int)
y_index = np.linspace(0, y_points.size-2, y_points.size-1, dtype=int)
z_index = np.linspace(0, z_points.size-2, z_points.size-1, dtype=int)

data_index = np.arange(np.size(x0))
matrix = np.concatenate([[-1, -1, -1], data_index])

for i_z in z_index:
    for i_y in y_index:
        for i_x in x_index:

            margin = x0*0 + initial_margin
            delta_margin = x0*0 + initial_delta_margin

            # print('index', i_x, i_y, i_z)
            line_length = np.zeros(x0.size)

            receiver_in_x_flag = np.logical_and(
                    x0 < x_points[i_x + 1],
                    x0 > x_points[i_x]
                    )
            receiver_in_y_flag = np.logical_and(
                    y0 < y_points[i_y + 1],
                    y0 > y_points[i_y]
                    )
            receiver_in_z_flag = np.logical_and(
                    z0 < z_points[i_z + 1],
                    z0 > z_points[i_z]
                    )
            receiver_in_voxel_flag = np.logical_and(
                    np.logical_and(receiver_in_x_flag, receiver_in_y_flag),
                    receiver_in_z_flag)

            max_permitted_intersections = (
                    2 - receiver_in_voxel_flag.astype(float))
            min_permitted_intersections = receiver_in_voxel_flag.astype(float)

            intersection_alpha = np.column_stack((
                x_alpha[:, i_x], x_alpha[:, i_x + 1],
                y_alpha[:, i_y], y_alpha[:, i_y + 1],
                z_alpha[:, i_z], z_alpha[:, i_z + 1],
                ))

            intersection_x = x0[:, None] + x1[:, None]*intersection_alpha
            intersection_y = y0[:, None] + y1[:, None]*intersection_alpha
            intersection_z = z0[:, None] + z1[:, None]*intersection_alpha

            flag_x = np.logical_and(
                    intersection_x <= x_points[i_x + 1] + margin[:, None],
                    intersection_x >= x_points[i_x] - margin[:, None]
                    )
            flag_y = np.logical_and(
                    intersection_y <= y_points[i_y + 1] + margin[:, None],
                    intersection_y >= y_points[i_y] - margin[:, None]
                    )
            flag_z = np.logical_and(
                    intersection_z <= z_points[i_z + 1] + margin[:, None],
                    intersection_z >= z_points[i_z] - margin[:, None]
                    )
            tmpflag = (
                    flag_x.astype(float) +
                    flag_y.astype(float) +
                    flag_z.astype(float))
            flag = tmpflag > 2

            negative_alpha_margin = margin/(sr*np.sin(elevation * np.pi/180))
            flag = np.logical_and(
                    flag,
                    intersection_alpha + negative_alpha_margin[:, None] > 0)

            min_permitted_intersections[
                    np.logical_and(
                        flag.sum(axis=1) > 0,
                        np.logical_not(receiver_in_voxel_flag)
                        )
                    ] = 2

            pre_iteration_flag = flag

            negated_flag = np.zeros(np.size(x0)).astype(bool)

            iteration_index_flag = np.logical_or(
                    flag.sum(axis=1) > max_permitted_intersections,
                    flag.sum(axis=1) < min_permitted_intersections
                    )
            while np.logical_or(
                    flag.sum(axis=1) > max_permitted_intersections,
                    flag.sum(axis=1) < min_permitted_intersections
                    ).sum() > 0:

                decrease_flag = flag.sum(axis=1) > max_permitted_intersections
                increase_flag = flag.sum(axis=1) < min_permitted_intersections

                # change
                delta_margin = delta_margin * (
                        1 + (delta_margin_factor - 1)*(np.logical_and(
                                decrease_flag, negated_flag
                            )).astype(float)
                        )

                margin = (
                        margin
                        - delta_margin*(decrease_flag).astype(float)
                        + delta_margin*(increase_flag).astype(float)
                        )

                negated_flag = increase_flag
                # print('margin', margin)

                flag_x = np.logical_and(
                        intersection_x < x_points[i_x + 1] + margin[:, None],
                        intersection_x > x_points[i_x] - margin[:, None]
                        )

                flag_y = np.logical_and(
                        intersection_y < y_points[i_y + 1] + margin[:, None],
                        intersection_y > y_points[i_y] - margin[:, None]
                        )

                flag_z = np.logical_and(
                        intersection_z < z_points[i_z + 1] + margin[:, None],
                        intersection_z > z_points[i_z] - margin[:, None]
                        )

                tmpflag = (
                        flag_x.astype(float) +
                        flag_y.astype(float) +
                        flag_z.astype(float))
                flag = tmpflag > 2

                negative_alpha_margin = margin/(
                        sr*np.sin(elevation * np.pi/180))
                flag = np.logical_and(
                        flag,
                        (intersection_alpha +
                            negative_alpha_margin[:, None]) > 0)

                # print(margin[iteration_index_flag],
                #       iteration_index_flag.sum())
                # print(pre_iteration_flag[iteration_index_flag])
                # print(flag[iteration_index_flag])

            # print(margin)

            # di -> double intersection
            # si -> single intersection
            # ni -> no intersection
            # i -> index
            # f -> flag
            # a -> alpha
            di_i = data_index[flag.sum(axis=1) == 2]
            si_i = data_index[flag.sum(axis=1) == 1]
            ni_i = data_index[flag.sum(axis=1) == 0]

            di_f = flag[di_i, :]
            si_f = flag[si_i, :]
            ni_f = flag[ni_i, :]

            di_a = intersection_alpha[di_i, :]
            si_a = intersection_alpha[si_i, :]
            ni_a = intersection_alpha[ni_i, :]

            ####
            di_a_straightened = np.reshape(di_a, di_a.size)
            di_f_straightened = np.reshape(di_f, di_f.size)
            di_a_selected = di_a_straightened[di_f_straightened]
            di_a_2d = np.reshape(di_a_selected,
                                 (int(di_a_selected.size / 2), 2))

            di_a_delta = di_a_2d[:, 1] - di_a_2d[:, 0]
            line_length[di_i] = np.sqrt(
                    (x1[di_i]*di_a_delta)**2 +
                    (y1[di_i]*di_a_delta)**2 +
                    (z1[di_i]*di_a_delta)**2)

            si_a_straightened = np.reshape(si_a, si_a.size)
            si_f_straightened = np.reshape(si_f, si_f.size)
            si_a_selected = si_a_straightened[si_f_straightened]

            line_length[si_i] = np.sqrt(
                    (si_a_selected*x1[si_i])**2 +
                    (si_a_selected*y1[si_i])**2 +
                    (si_a_selected*z1[si_i])**2)

            line_length = np.concatenate([[i_x, i_y, i_z], line_length])
            matrix = np.column_stack((matrix, line_length))

            post_iteration_flag = flag

            # diagnostics that will be printed for inspection
            diff_flag = np.not_equal(pre_iteration_flag, post_iteration_flag)
            diff_index = data_index[diff_flag.sum(axis=1) != 0]
            print(np.not_equal(pre_iteration_flag, post_iteration_flag).sum(),
                  diff_index,
                  pre_iteration_flag[diff_index, :],
                  post_iteration_flag[diff_index, :])

            # if plotting is required
            if(plot_epoch_index >= 0 and
                    min_permitted_intersections[plot_epoch_index] > 0):

                for i in [0, 1]:
                    tmpy = np.linspace(y_points[i_y], y_points[i_y+1], 2)
                    tmpz = np.linspace(z_points[i_z], z_points[i_z+1], 2)
                    tmpY, tmpZ = np.meshgrid(tmpy, tmpz)
                    tmpX = tmpY*0 + x_points[i_x + i]
                    ax.plot_wireframe(tmpX, tmpY, tmpZ)
                for i in [0, 1]:
                    tmpx = np.linspace(x_points[i_x], x_points[i_x+1], 2)
                    tmpz = np.linspace(z_points[i_z], z_points[i_z+1], 2)
                    tmpX, tmpZ = np.meshgrid(tmpx, tmpz)
                    tmpY = tmpX*0 + y_points[i_y + i]
                    ax.plot_wireframe(tmpX, tmpY, tmpZ)
                for i in [0, 1]:
                    tmpx = np.linspace(x_points[i_x], x_points[i_x+1], 2)
                    tmpy = np.linspace(y_points[i_y], y_points[i_y+1], 2)
                    tmpX, tmpY = np.meshgrid(tmpx, tmpy)
                    tmpZ = tmpX*0 + z_points[i_z + i]
                    ax.plot_wireframe(tmpX, tmpY, tmpZ)

                tmpalpha = np.linspace(
                        intersection_alpha[
                            plot_epoch_index, flag[plot_epoch_index, :]
                            ][0] * (
                            ~receiver_in_voxel_flag[plot_epoch_index]
                            ).astype(float),

                        intersection_alpha[
                            plot_epoch_index, flag[plot_epoch_index, :]
                            ][-1],

                        2)
                tmpx = x0[plot_epoch_index] + tmpalpha*x1[plot_epoch_index]
                tmpy = y0[plot_epoch_index] + tmpalpha*y1[plot_epoch_index]
                tmpz = z0[plot_epoch_index] + tmpalpha*z1[plot_epoch_index]
                ax.scatter(tmpx, tmpy, tmpz, c='r')

                if(min_permitted_intersections[plot_epoch_index] == 1):
                    tmppei = plot_epoch_index
                    tmpalpha1 = intersection_alpha[tmppei, flag[tmppei, :]][0]
                    print('BEGIN')
                    print(i_x, i_y, i_z)
                    print(x0[tmppei] + tmpalpha1*x1[tmppei])
                    print(y0[tmppei] + tmpalpha1*y1[tmppei])
                    print(z0[tmppei] + tmpalpha1*z1[tmppei])
                    print('END')
                if(min_permitted_intersections[plot_epoch_index] == 2):
                    tmppei = plot_epoch_index
                    tmpalpha = intersection_alpha[tmppei, flag[tmppei, :]]
                    print('BEGIN')
                    print(i_x, i_y, i_z)
                    print(x0[tmppei] + tmpalpha.min()*x1[tmppei])
                    print(y0[tmppei] + tmpalpha.min()*y1[tmppei])
                    print(z0[tmppei] + tmpalpha.min()*z1[tmppei])
                    print(x0[tmppei] + tmpalpha.max()*x1[tmppei])
                    print(y0[tmppei] + tmpalpha.max()*y1[tmppei])
                    print(z0[tmppei] + tmpalpha.max()*z1[tmppei])
                    print('END')

                plt.show(block=False)


np.savetxt(output_matrix_file_name, matrix)

# --------------------------------------------------------
# | for doing the validation calculation with whole space as one voxel
# --------------------------------------------------------
validation_x_alpha = np.column_stack([x_alpha[:, 0], x_alpha[:, -1]])
validation_y_alpha = np.column_stack([y_alpha[:, 0], y_alpha[:, -1]])
validation_z_alpha = np.column_stack([z_alpha[:, 0], z_alpha[:, -1]])
validation_x_points = np.array([x_points[0], x_points[-1]])
validation_y_points = np.array([y_points[0], y_points[-1]])
validation_z_points = np.array([z_points[0], z_points[-1]])


data_index = np.arange(np.size(x0))

margin = x0*0 + initial_margin
delta_margin = x0*0 + initial_delta_margin

line_length = np.zeros(x0.size)
validation_x1 = line_length*0
validation_y1 = line_length*0
validation_z1 = line_length*0
validation_x2 = line_length*0
validation_y2 = line_length*0
validation_z2 = line_length*0
validation_no_of_int = line_length*0

receiver_in_x_flag = np.logical_and(
        x0 < validation_x_points[1],
        x0 > validation_x_points[0]
        )
receiver_in_y_flag = np.logical_and(
        y0 < validation_y_points[1],
        y0 > validation_y_points[0]
        )
receiver_in_z_flag = np.logical_and(
        z0 < validation_z_points[1],
        z0 > validation_z_points[0]
        )
receiver_in_voxel_flag = np.logical_and(
        np.logical_and(receiver_in_x_flag, receiver_in_y_flag),
        receiver_in_z_flag)

max_permitted_intersections = 2 - receiver_in_voxel_flag.astype(float)
min_permitted_intersections = receiver_in_voxel_flag.astype(float)

intersection_alpha = np.column_stack((
    validation_x_alpha[:, 0], validation_x_alpha[:, 1],
    validation_y_alpha[:, 0], validation_y_alpha[:, 1],
    validation_z_alpha[:, 0], validation_z_alpha[:, 1],
    ))

intersection_x = x0[:, None] + x1[:, None]*intersection_alpha
intersection_y = y0[:, None] + y1[:, None]*intersection_alpha
intersection_z = z0[:, None] + z1[:, None]*intersection_alpha

flag_x = np.logical_and(
        intersection_x <= validation_x_points[1] + margin[:, None],
        intersection_x >= validation_x_points[0] - margin[:, None]
        )
flag_y = np.logical_and(
        intersection_y <= validation_y_points[1] + margin[:, None],
        intersection_y >= validation_y_points[0] - margin[:, None]
        )
flag_z = np.logical_and(
        intersection_z <= validation_z_points[1] + margin[:, None],
        intersection_z >= validation_z_points[0] - margin[:, None]
        )
tmpflag = flag_x.astype(float) + flag_y.astype(float) + flag_z.astype(float)
flag = tmpflag > 2

negative_alpha_margin = margin/(sr*np.sin(elevation * np.pi/180))
flag = np.logical_and(
        flag,
        (intersection_alpha + negative_alpha_margin[:, None]) > 0)

min_permitted_intersections[
        np.logical_and(
            flag.sum(axis=1) > 0,
            np.logical_not(receiver_in_voxel_flag)
            )
        ] = 2


pre_iteration_flag = flag

negated_flag = np.zeros(np.size(x0)).astype(bool)

iteration_index_flag = np.logical_or(
        flag.sum(axis=1) > max_permitted_intersections,
        flag.sum(axis=1) < min_permitted_intersections
        )
while np.logical_or(
        flag.sum(axis=1) > max_permitted_intersections,
        flag.sum(axis=1) < min_permitted_intersections
        ).sum() > 0:

    decrease_flag = flag.sum(axis=1) > max_permitted_intersections
    increase_flag = flag.sum(axis=1) < min_permitted_intersections

    delta_margin = delta_margin*(
            1 + (delta_margin_factor - 1)*(np.logical_and(
                    decrease_flag, negated_flag
                )).astype(float)
            )

    margin = (
            margin
            - delta_margin*(decrease_flag).astype(float)
            + delta_margin*(increase_flag).astype(float)
            )

    negated_flag = increase_flag

    flag_x = np.logical_and(
            intersection_x < validation_x_points[1] + margin[:, None],
            intersection_x > validation_x_points[0] - margin[:, None]
            )
    flag_y = np.logical_and(
            intersection_y < validation_y_points[1] + margin[:, None],
            intersection_y > validation_y_points[0] - margin[:, None]
            )
    flag_z = np.logical_and(
            intersection_z < validation_z_points[1] + margin[:, None],
            intersection_z > validation_z_points[0] - margin[:, None]

            )

    tmpflag = (
            flag_x.astype(float) +
            flag_y.astype(float) +
            flag_z.astype(float))
    flag = tmpflag > 2

    negative_alpha_margin = margin/(sr*np.sin(elevation * np.pi/180))
    flag = np.logical_and(
            flag,
            (intersection_alpha + negative_alpha_margin[:, None]) > 0)


post_iteration_flag = flag

# diagnostics of validation calculation for inspection
diff_flag = np.not_equal(pre_iteration_flag, post_iteration_flag)
diff_index = data_index[diff_flag.sum(axis=1) != 0]
print('validate:', np.not_equal(pre_iteration_flag, post_iteration_flag).sum(),
      diff_index,
      pre_iteration_flag[diff_index, :],
      post_iteration_flag[diff_index, :])

# di -> double intersection
# si -> single intersection
# ni -> no intersection
# i -> index
# f -> flag
# a -> alpha
di_i = data_index[flag.sum(axis=1) == 2]
si_i = data_index[flag.sum(axis=1) == 1]
ni_i = data_index[flag.sum(axis=1) == 0]

di_f = flag[di_i, :]
si_f = flag[si_i, :]
ni_f = flag[ni_i, :]

di_a = intersection_alpha[di_i, :]
si_a = intersection_alpha[si_i, :]
ni_a = intersection_alpha[ni_i, :]

####
di_a_straightened = np.reshape(di_a, di_a.size)
di_f_straightened = np.reshape(di_f, di_f.size)
di_a_selected = di_a_straightened[di_f_straightened]
di_a_2d = np.reshape(di_a_selected, (int(di_a_selected.size / 2), 2))

di_a_delta = di_a_2d[:, 1] - di_a_2d[:, 0]
line_length[di_i] = np.sqrt(
        (x1[di_i]*di_a_delta)**2 +
        (y1[di_i]*di_a_delta)**2 +
        (z1[di_i]*di_a_delta)**2)
validation_x1[di_i] = x0[di_i] + x1[di_i]*di_a_2d[:, 0]
validation_y1[di_i] = y0[di_i] + y1[di_i]*di_a_2d[:, 0]
validation_z1[di_i] = z0[di_i] + z1[di_i]*di_a_2d[:, 0]
validation_x2[di_i] = x0[di_i] + x1[di_i]*di_a_2d[:, 1]
validation_y2[di_i] = y0[di_i] + y1[di_i]*di_a_2d[:, 1]
validation_z2[di_i] = z0[di_i] + z1[di_i]*di_a_2d[:, 1]

si_a_straightened = np.reshape(si_a, si_a.size)
si_f_straightened = np.reshape(si_f, si_f.size)
si_a_selected = si_a_straightened[si_f_straightened]

line_length[si_i] = np.sqrt(
        (si_a_selected*x1[si_i])**2 +
        (si_a_selected*y1[si_i])**2 +
        (si_a_selected*z1[si_i])**2)
validation_x1[si_i] = x0[si_i]
validation_y1[si_i] = y0[si_i]
validation_z1[si_i] = z0[si_i]
validation_x2[si_i] = x0[si_i] + x1[si_i]*si_a_selected
validation_y2[si_i] = y0[si_i] + y1[si_i]*si_a_selected
validation_z2[si_i] = z0[si_i] + z1[si_i]*si_a_selected

validation_z1[ni_i] = -1
validation_z2[ni_i] = -1

# | This helps us understand weather GPS receiver is inside the grid space or
# | outside. If outside then weather rays are intersecting with grid space or
# | not.
validation_no_of_int[si_i] = 1
validation_no_of_int[di_i] = 2
validation_no_of_int[ni_i] = 0

validation_data = np.column_stack([data_index, line_length])  # 0 and 1
validation_data = np.column_stack([validation_data, validation_no_of_int])  # 2
validation_data = np.column_stack([validation_data, validation_x1])  # 3
validation_data = np.column_stack([validation_data, validation_y1])  # 4
validation_data = np.column_stack([validation_data, validation_z1])  # 5
validation_data = np.column_stack([validation_data, validation_x2])  # 6
validation_data = np.column_stack([validation_data, validation_y2])  # 7
validation_data = np.column_stack([validation_data, validation_z2])  # 8

np.savetxt(validation_file_name, validation_data)


# --------------------------------------------------------
# | for some plotting
# --------------------------------------------------------
if master_plot_flag:

    unique_index = np.unique(receiver_number, return_index=True)[1]
    mpluse('fivethirtyeight')
    tmpfig = plt.figure('map like')
    ax = plt.subplot(111)
    tmplon = x_points/(x_deg_to_meter_factor) + min_lon
    tmplat = y_points/(y_deg_to_meter_factor) + min_lat
    tmpLon, tmpLat = np.meshgrid(tmplon, tmplat)
    ax.plot(tmpLon, tmpLat, '#0F95D7')
    tmpLat, tmpLon = np.meshgrid(tmplat, tmplon)
    ax.plot(tmpLon, tmpLat, '#0F95D7')

    ax.scatter(receiver_lon[unique_index], receiver_lat[unique_index],
               c='#FF2700', s=100)
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    ax.set_title('Voxel grid and Receivers')
    # tmpfig.tight_layout()

    tmpfig_grid = plt.figure('grid')
    ax = plt.axes(projection=Axes3D.name)
    ax.set_xlabel('Along longitude (m)', labelpad=20)
    ax.set_ylabel('Along latitude (m)', labelpad=20)
    ax.set_zlabel('Geodetic height (m)', labelpad=20)
    ax.set_title('Voxel space')

    tmpx = x_points
    tmpy = y_points
    tmpX, tmpY = np.meshgrid(tmpx, tmpy)
    tmpZ = tmpX*0 + z_points.min()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)
    tmpx = np.linspace(x_points.min(), x_points.max(), 2)
    tmpy = np.linspace(y_points.min(), y_points.max(), 2)
    tmpX, tmpY = np.meshgrid(tmpx, tmpy)
    tmpZ = tmpX*0 + z_points.max()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)

    tmpx = x_points
    tmpz = z_points
    tmpX, tmpZ = np.meshgrid(tmpx, tmpz)
    tmpY = tmpX*0 + y_points.max()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)
    tmpx = np.linspace(x_points.min(), x_points.max(), 2)
    tmpz = np.linspace(z_points.min(), z_points.max(), 2)
    tmpX, tmpZ = np.meshgrid(tmpx, tmpz)
    tmpY = tmpX*0 + y_points.min()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)

    ax.scatter(x0[unique_index], y0[unique_index], z0[unique_index],
               c='r', s=80, label='Receiver locations')
    ax.legend()

    # mpluse('default')

    tmpfig = plt.figure('rays')
    ax = plt.axes(projection=Axes3D.name)
    ax.set_xlabel('Longitude (m)', labelpad=20)
    ax.set_ylabel('Latitude (m)', labelpad=20)
    ax.set_zlabel('Geodetic height (m)', labelpad=20)
    ax.set_title('Rays , Ray Cutoff Height =' +
                 str(ray_cutoff_height/1000) + 'km')

    tmpx = x_points
    tmpy = y_points
    tmpX, tmpY = np.meshgrid(tmpx, tmpy)
    tmpZ = tmpX*0 + z_points.min()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)
    tmpx = np.linspace(x_points.min(), x_points.max(), 2)
    tmpy = np.linspace(y_points.min(), y_points.max(), 2)
    tmpX, tmpY = np.meshgrid(tmpx, tmpy)
    tmpZ = tmpX*0 + z_points.max()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)

    tmpx = x_points
    tmpz = z_points
    tmpX, tmpZ = np.meshgrid(tmpx, tmpz)
    tmpY = tmpX*0 + y_points.max()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)
    tmpx = np.linspace(x_points.min(), x_points.max(), 2)
    tmpz = np.linspace(z_points.min(), z_points.max(), 2)
    tmpX, tmpZ = np.meshgrid(tmpx, tmpz)
    tmpY = tmpX*0 + y_points.min()
    ax.plot_wireframe(tmpX, tmpY, tmpZ)

    ax.scatter(x0[unique_index], y0[unique_index], z0[unique_index],
               c='r', s=80, label='Receiver locations')
    ax.legend()

    # | plotting rays of each receiver individually
    for r_i in np.unique(receiver_number):
        tmpflag = receiver_number == r_i
        tmpx = np.vstack((
            validation_x1[tmpflag],
            validation_x2[tmpflag],
            validation_x1[tmpflag])
            ).reshape((-1), order='F')

        tmpy = np.vstack((
            validation_y1[tmpflag],
            validation_y2[tmpflag],
            validation_y1[tmpflag])
            ).reshape((-1), order='F')

        tmpz = np.vstack((
            validation_z1[tmpflag],
            validation_z2[tmpflag],
            validation_z1[tmpflag])
            ).reshape((-1), order='F')

        ax.plot(tmpx, tmpy, tmpz, linewidth=.5, color='#ff7f0e')

    plt.show()
