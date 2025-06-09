# Swarming Toolpath Generator

The main code is

    particle_traj_opt_casadi_QP_bounds.py

It will require a few common packages (see beginning of code). In particular CASADI is needed for optimization.

### Inputs
The current code runs with an example on an open hole tensile specimen.

The code requires an FEM simulation of a load case from a (3 x 3 x n) array (currently in 'FEM_results.mat') and the corresponding mesh (currently in 'std_specimen.STL')

If changing of part, the starting locations of agents should be adapted. They are currently given as

    start_grid = np.transpose(np.meshgrid(np.arange(0.5, 36.5, desired_distance), [1], [61])).reshape(-1, 3)


### 2d Implementation
The current version only runs in 2D cases. However, the FEM and mesh are 3D.
The code will consider the stress obtained on the plane Y=1 and project it on this X,Z plane. It allows particles only to move in the X,Z plane with Y=1.

Since the plan was to extend it to 3D, all code lines that are a 2D adaptation are marked with a comment saying '2D version'. They can easily be found by searchin in the text.

### Algorithm details
All mathematical details of functions can be found in the paper.

The main tuning handle is

    K_opt = 5

Lower values (e.g. K_opt = 0.5) produce more uniform lines, while higher values (e.g. K_opt = 50) produce better stress alignment.

Line spacing can be set width

    desired_distance = 0.4

If particles don't reach the opposite end of the part, the number of steps should simply be increased

    n_steps = 100

## Other functions

The following functions can be activated or deactivated at the bottom of the code by toggling True/False condtions
- plotting
- video creation
- stress alignment score
- computing and saving line distances (which can be plotted by a dedicated code, see open hole results folder)
