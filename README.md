# tropotomo

GNSS receiver network tropospheric tomography matrix calculator. For details regarding detailed explanation read `explanation.pdf`. Calculations done for small scale GNSS network with perfectly cuboidal voxels.

Basic tomography matrix equation is **y=Ax**, where **y** is a column matrix containing all the line integrated measurements, **x** is a column matrix of unknown parameters in each voxel and **A** is a rectangular matrix containing the intersection lengths of each line/ray and voxel. Each row in **A** corresponds to a ray and each column corresponds to voxel.

The goal of `intersection_matrix.py` is to calculate the matrix **A**.

## Usage

`intersection_matrix.py` does the ray and voxel intersection length calculations. Before running this script, `tomo_config.py` should be edited accordingly.

`tomo_config.py` will have different variables that will be imported in `intersection_matrix.py`. For integrating `intersection_matrix.py` into different scripts it is enough to make sure all the required variables are defined and imported.

To run use the command `python3 intersection_matrix.py`.

### Output messages

When running the code, it will be printing some numbers. Each line is for one particular voxel. There is no loop exclusively over each and every ray. Each ray corresponds to a row in data and all the variables. Most calculations are written such that there is no need to loop over each ray by array operations even if they get complicated.

There will be 6 flags for each combination of ray and voxel intersection. These 6 flags are for 6 (surface) planes of voxel. Ideally 2 of these flags should be True if ray is actually intersecting the voxel. If there is no intersection then all of them will be False. If the GNSS receiver is inside the voxel then any one of the flag will always be True. Once the correct flags are set, corresponding intersection points can be used to find intersection length.

During the calculation the number of flags may not be exactly right in the first time for some rays. In such cases the correct flags will be obtained after some number of iterations. These changes will be printed by code.

There will be 4 columns, printed for each voxel in a line.
* 1 - an integer, the total number of flags changed for 6 × number of rays.
* 2 - a 1D array of indexes (starting from 0) of rays whose flags were changed.
* 3 - a 2D array of pre iteration flags of each ray that had changes in each row and 6 flags in each column.
* 4 - a 2D array of post iteration flags of each ray that had changes in each row and 6 flags in each column.

Arrays are identified by matching '[' and ']' and order of flags is same as array of indexes.
Order of 6 flags is intersection of ray with surfaces of min\_x, max\_x, min\_y, max\_y, min\_z, max\_z of the voxel.

### Data files and format

#### Input file

All data files are of plane text. Input data file should have columns of receiver numbers, receiver longitude, latitude, altitude (from reference ellipsoid), azimuth and elevation angle. Column number of these data should be given in `tomo_config.py`. Any number of columns can be present in input data file, but any column other than these will be ignored. Number of lines in the beginning that should not be considered or number of header lines that should not be considered should be given.

#### Output files

Two files will be output by `intersection_matrix.py`. One is the actual matrix **A** and other is a validation data file. Validation data file has the same calculations at **A** but with the whole voxel space as a single voxel. Sum of individual lengths of any ray should be equal to the length calculated in validation data file. This can be used to inspect if there are any errors in calculations.

**A** matrix file has 3 rows in the beginning and one column at the beginning which give info regarding the ray and voxel combination. The three rows give the x, y and z voxel indices (from 0) of the voxel, first column will give the index of the ray. Values in the top left corner will be filled with -1.

Validation data file will have ray index, intersection length, number of intersection points, starting intersection point location x, y, z from receiver and ending intersection point location x, y, z towards the satellite. x, y, z values here are with respect the axes defined along the edges of voxel space with origin at corner. These intersection points will be on the surfaces of the whole voxel space, and starting points will be the receiver location if receiver is inside the voxel space.

Note that the total length from **A** and validation data may not be exactly equal due to loss of precision while saving data in output files.

To change the output file format edit the `np.savetxt(file_name_str, numpy_nD_array)` line in `intersection_matrix.py` file.

## sample data generator
The script `make_sample_data.py` can be used to make a sample input data.
