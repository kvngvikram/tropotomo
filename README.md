# tropotomo
GNSS receiver network tropospheric tomographic matrix calculator.


When running the code, it will be printeng some numbers. Each line is for one particular voxel. There is no loop exclusively over each and every ray. Each ray corresponds to a row in data and all the variables. Most calculations are written such that there is no need to loop over each ray by array operations even if they get complicated.

There will be 6 flags for each combination of ray and voxel intersection. These 6 flags are for 6 (surface) planes of voxel. Ideally 2 of these flags should be True if ray is actually intersecting the voxel. If there is no intersection then all of them will be False. If the GNSS receiver is inside the voxel then any one of the flag will always be True. Once the correct flags are set then corresponding intersection points can be used to find intersection length.

During the calulation the number of flags may not be exactly right in the first time for some rays. In such cases the correct flags will be obtained after some number of iterations. These changes will be printed by code.

There will be 4 columns, printed for each voxel in a line.
* 1 - an integer, the total number of flags changed for 6 * number of rays.
* 2 - a 1D array of indexes (starting from 0) of rays whose flags were changed.
* 3 - a 2D array of pre iteration flags of each ray that had changes in each row and 6 flags in each column.
* 4 - a 2D array of post iteration flags of each ray that had changes in each row and 6 flags in each column.

Arrays are identified by matching '[' and ']' and order of flags is same as array of indexes.
Order of 6 flags is intersection with surface of min_x, max_x, min_y, max_y, min_z, max_z.

## sample data generator
A seperate code is written just to generate a sample data that can be used for testing.
