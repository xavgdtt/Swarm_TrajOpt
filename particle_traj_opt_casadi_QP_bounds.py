# -*- coding: utf-8 -*-
"""
@author: Xavier Guidetti

"""

%matplotlib auto

import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from matplotlib.cm import turbo
from mpl_toolkits.mplot3d import Axes3D
import time
from casadi import *
import trimesh
import cv2
import glob

"""
Functions definitions
"""

# ----------------------
# 1. Functions for stress
# ----------------------

def interpolateStress(query_points,interpolator):
    """
    Interpolate stress at specified spatial locations by nearest node method.

    """
        
    IntrpStress_swap = interpolator(query_points)
    IntrpStress = IntrpStress_swap.swapaxes(0,2)
        
    return IntrpStress

def computeMaxPrincipalStress(stress_tensors):
    
       
    """
    Compute maximum principal stress vector and corresponding stress magnitude from stress tensors.

    Args:
        stress_tensors (array_like): Array of stress tensors with shape (3, 3, n), where n is the number of stress tensors.

    Returns:
        PSmax (array): Array of maximum principal stress vectors with shape (n, 3), where n is the number of stress tensors.
        SigmaMax (array): Array of corresponding stress magnitude with shape (n), where n is the number of stress tensors.
    """

    n = stress_tensors.shape[2]  # Number of stress tensors
    PSmax = np.zeros((n, 3))  # Array to store principal stress vectors
    SigmaMax = np.zeros((n))  # Array to store sorted eigenvalues

    for i in range(n):
        cst = stress_tensors[:, :, i]  # Extract stress tensor for the i-th stress tensor
        D, V = np.linalg.eig(cst)  # Compute eigenvalues and eigenvectors of the stress tensor
        idx = np.argmax(np.abs(D))  # Find largest eigenvalue by absolute value
        SigmaMax[i] = D[idx]  # Store the largest eigenvalue
        PSmax[i, :] = V[:, idx].T  # Store the eigenvector corresponding to the maximum eigenvalue

    return PSmax, SigmaMax

def align_psmax(PSmax):
    """
    Aligns the rows of PSmax matrix in a way that the third element of each row is non-negative.

    Args:
        PSmax (numpy array): A 2D numpy array representing the PSmax matrix.

    Returns:
        numpy array: The aligned PSmax matrix.
        
 
    """
    
    is_not_aligned = PSmax[:, 2] < 0
    PSmax[is_not_aligned, :] = -PSmax[is_not_aligned, :]
    return PSmax


# ----------------------
# 2. Functions for distance computation
# ----------------------

def compute_orth_dist(points, points_old):
    """
    Calculates the orthogonal distance between points in a 2D space.
    Uses the previous motion of agents as a reference direction to obtain the
    reference axis along which the orthogonal direction is computed
    
    Parameters:
        - points: ndarray, shape (n, 3)
            Array containing the current 3D points.
        - points_old: ndarray, shape (n, 3)
            Array containing the previous 3D points.

    """
    
    d = np.zeros(points.shape[0]-1) # Initialize array to store distance vectors

    for i in range(points.shape[0]-1): # iterate over all agents
        j = i + 1 # compute distance between agent and following agent

        if np.isnan(points_old[i,:]).all() and np.isnan(points_old[j,:]).all():
            d[i] = np.linalg.norm( points[i,:] - points[j,:] )
        else:
            if np.isnan(points_old[i,:]).all():
                va = points[j,:] - points_old[j,:]
                vb = points[i,:] - points_old[j,:]
                ve = points[i,:] - points[j,:]
            else:
                va = points[i,:] - points_old[i,:]
                vb = points[j,:] - points_old[i,:]
                ve = points[j,:] - points[i,:]

            d[i] = np.linalg.norm( vb - (dot(va,vb)*va / dot(va,va)) ) # 2D version : check if this applies in 3D, it might not
           
    return d


# ----------------------
# 3. Functions for energies
# ----------------------

def compute_alignment(points, points_old, x_des):
    """
    Calculates the alignment score between agents
    
    Details of alignment formula are in the paper
    
    """
        
    al = MX.zeros(points.shape[0]-1) # Initialize array to store distances

    for i in range(points.shape[0]-1): # iterate over all agents
        j = i + 1 # compute alignment between agent and following agent
        
        if np.isnan(points_old[i,:,:]).any() and np.isnan(points_old[j,:,:]).any(): 
            
            if np.isnan(points_old[i,:,1]).any():
                dmean = x_des[j,:] - points_old[j,:,1]
            else:
                dmean = x_des[i,:] - points_old[i,:,1]
                         
            #The case where both particles have no past is not covered yet!
            
            v = points[j,:] - points[i,:]
            
            d = [dmean[1], -dmean[0]] # 2D version : must computed perpendiculars differently in 3D
                        
            
            d = d/np.linalg.norm(d)
            d = np.expand_dims(d,axis=0)
            proj = (dot(v,d) - desired_distance)**2
            oproj = v - dot(v,d)*d
            rej = dot(oproj,oproj)
            
            al[i] = proj + rej
                    
        else:
            # this is the case where both particles have a past

            v = points[j,:] - points[i,:]
            
            if np.isnan(points_old[i,:,:]).any():
                dmean = points_old[j,:,1] - points_old[j,:,0]
                
            else: 
                if np.isnan(points_old[j,:,:]).any():
                    dmean = points_old[i,:,1] - points_old[i,:,0]
                    
                else:
                    # formulation 1
                    di = points_old[i,:,1] - points_old[i,:,0]
                    dj = points_old[j,:,1] - points_old[j,:,0]
                    dmean = (di+dj)/2
                            
            d = [dmean[1], -dmean[0]] # 2D version : must make perpendicular differently in 3D
                 
          
            d = d/np.linalg.norm(d)
            d = np.expand_dims(d,axis=0)
            proj = (dot(v,d) - desired_distance)**2
            oproj = v - dot(v,d)*d
            rej = dot(oproj,oproj)
            
            al[i] = proj + rej       
           
    return al


def potential(x, x_old, x_des):
    """
    Calculates the spacing potential for the entire set of points.

    Returns:
        - J: float
            The calculated spacing potential.
    """

    # Compute all alignments between points
    indiv_alignments = compute_alignment(x,x_old,x_des)

    # Merges all individual alignments into one score
    J = sum1(indiv_alignments)
    
    return J

def control_cost(x, x_des, x_old, m):
    """
    Calculates the control cost for the particles optimization problem.

    Details of formula in the paper
    """
    
    dist = x_des - x
    square_of_dist = (sum2(dist**2))
    
    ctrl_cost = sum1(m*(square_of_dist))
    
    return ctrl_cost


def get_location_from_dev(x_dev, x_des, x_old, step_size, opt_loop):
    """
    Computes a point absolute location based on the local deviation from the point projected 
    from the past motion.
    If past motion is absent, it does so based on deviation from desired point
    from stress field
    
    -step_size <= x_dev <= step_size
    """
    
    ## if they have no past, they must be processed
    past_missing = np.any(np.any(np.isnan(x_old),axis=2),axis=1)
    
    points_number = x_old.shape[0]
    
    R = np.zeros([points_number,2,2])
    x_momentum = np.zeros([points_number,2])
    
    if past_missing.sum() == points_number: # first steps where all elements have no past    
        # just make them all go straight
        R[past_missing,:,:] = [[1,0],[0,1]]
        x_momentum[past_missing,:] = x_des[past_missing,:] #reference at desired
    else: #if neighbots have past, interpolate the past positions of the the neighbors
        indices = np.arange(points_number)
        x_old[:,0,0] =  np.interp(indices, indices[~past_missing], x_old[~past_missing,0,0])
        x_old[:,1,0] =  np.interp(indices, indices[~past_missing], x_old[~past_missing,1,0])
        x_old[:,0,1] =  np.interp(indices, indices[~past_missing], x_old[~past_missing,0,1])
        x_old[:,1,1] =  np.interp(indices, indices[~past_missing], x_old[~past_missing,1,1])
        
        previous_step = x_old[:,:,1] - x_old[:,:,0]
    
        ang = np.arctan2(previous_step[:,1],previous_step[:,0])
        R[:,0,0] = np.cos(ang)
        R[:,0,1] = -np.sin(ang)
        R[:,1,0] = np.sin(ang)
        R[:,1,1] = np.cos(ang) #2D VERSION : fix rotation matrices to change of reference system
        
        step_norm = np.linalg.norm(previous_step,axis=1)
        x_momentum = x_old[:,:,1] + previous_step/step_norm[:,np.newaxis] * step_size
    
    # different syntax depending on variable type
    if opt_loop:
        # shift from x_momentum
        x = MX.zeros(x_dev.shape)
        rot_dev = MX.zeros(x_dev.shape)
        # Rotation
        rot_dev[:,0] = sum2(x_dev * R[:,:,0])
        rot_dev[:,1] = sum2(x_dev * R[:,:,1])
    else:
        # shift from x_momentum
        x = np.zeros(x_dev.shape)
        rot_dev = np.zeros(x_dev.shape)
        # Rotation
        rot_dev[:,0] = np.sum(x_dev * R[:,:,0],axis=1)
        rot_dev[:,1] = np.sum(x_dev * R[:,:,1],axis=1)
    # Sum
    x = x_momentum + rot_dev
        
    return x

def tot_cost(x_dev, x_des, x_old, m, K, step_size):
    """
    Calculates the total cost function for a the particles optimization problem.

    Parameters:
        - x: ndarray, shape (n,)
            Array containing the current state variables.
        - x_des: ndarray, shape (n,)
            Array containing the desired state variables.
        - x_old: ndarray, shape (n,)
            Array containing the previous state variables.
        - m: xxx
            Array representing a weight parameter.
        - K: float
            Scalar value representing a weight parameter.

    Returns:
        - S: float
            Scalar value representing the total cost.
    """

    x = get_location_from_dev(x_dev, x_des, x_old, step_size, True)
    #x = x_dev + x_des

    # Calculate the control cost
    control = K * control_cost(x, x_des, x_old, m)

    # Calculate the spacing potential
    spacing = potential(x, x_old, x_des)
    
    #repulsion = ...

    # Calculate the total cost
    S = control + spacing # here a repulsion energy from specific locations can be added as an extra term + repulsion
    
    return S



# ----------------------
# 4. Functions for optimization
# ----------------------

def optimize_particles(trajs, mass, K_opt, j):
    """
    Optimizes the porticles positions for a given time step

    Parameters:
        - trajs: ndarray, shape (n, 3, m)
            Array containing the 3D trajectories for all agents.
        - mass: ndarray, shape (n,)
            Array containing the mass of all agents.
        - K_opt: float
            The optimization parmeter.
        - j: int
            The time step to optimize.

    Returns:
        - trajs: ndarray, shape (n, 3, m)
            The optimized 3D trajectories.
        - cost_per_part: cost of the solution
    """

    for_opt = np.all(~np.isnan(trajs[:, [0, 2], j+1]), axis=1) # exclude NaNs
    init_trajs = trajs[for_opt][:, [0, 2], j+1] # 2D version : constrain to planar
    desired_trajs = trajs[for_opt][:, [0, 2], j+1] # 2D version : constrain to planar
    if j == 0:
        previous_trajs = np.stack((trajs[for_opt][:, [0, 2], -1],trajs[for_opt][:, [0, 2], 0]),axis=2)
    else:
        previous_trajs = trajs[for_opt][:, [0, 2], j-1:j+1] # 2D version : constrain to planar
            
    init_dev = np.zeros(desired_trajs.shape)
    min_dev = -0.5*step_size
    max_dev = 0.5*step_size
    bounds_seq = (((min_dev,max_dev), ) * init_dev.shape[0])
    
    opti = casadi.Opti('conic')
    opts = {'osqp.verbose': 0}
    opti.solver('osqp',opts)
    #opti.solver('qpoases')
    

    x = opti.variable(desired_trajs.shape[0],desired_trajs.shape[1]) # deviation from going one step straight with momentum
    opti.minimize(tot_cost(x, desired_trajs, previous_trajs, mass[for_opt], K_opt, step_size))
    opti.subject_to( opti.bounded(1*min_dev,x[:,0],1*max_dev) ) # different bounds on different directions
    opti.subject_to( opti.bounded(0.5*min_dev,x[:,1],0.5*max_dev) ) # different bounds on different directions
    
    start_time = time.time() # time it
    sol = opti.solve()
    opt_duration.append(time.time() - start_time)
    
    optimizer = get_location_from_dev(sol.value(x), desired_trajs, previous_trajs, step_size, False)
    trajs[for_opt, :, j+1] = np.column_stack((optimizer[:,0], np.ones(optimizer.shape[0]), optimizer[:,1])) # 2D version : constrain to planar
    cost_per_part = sol.value(tot_cost(x, desired_trajs, previous_trajs, mass[for_opt], K_opt, step_size)) / for_opt.sum() # objective value divided by number of particles in the optimization

    return trajs, cost_per_part


# ----------------------
# 4. Functions for changing particles number
# ----------------------


def add_particle(trajs,mass,j,part_mesh):
    """
    Adds particle at the widest spacing
    """
    
    for_dist = np.all(~np.isnan(trajs[:, [0, 2], j + 1]), axis=1) # which particles to evaluate distances upon
    trajs_not_nan = trajs[for_dist, :, j + 1]
    trajs_not_nan_old = trajs[for_dist, :, j]
    point_dists_not_nan = compute_orth_dist(trajs_not_nan, trajs_not_nan_old)
    point_dists = np.zeros(trajs[:, :, j + 1].shape[0])
    have_diff = np.where(for_dist)[0]
    point_dists[have_diff[:-1]] = point_dists_not_nan
    spawn_ix_not_nan_all = np.argsort(point_dists_not_nan)[::-1]
    spawn_ix_all = np.argsort(point_dists)[::-1]

    for spawn_iter in range(point_dists_not_nan.shape[0]):
        spawn_ix = spawn_ix_all[spawn_iter]
        spawn_ix_not_nan = spawn_ix_not_nan_all[spawn_iter]
        potential_spawn = (trajs_not_nan[spawn_ix_not_nan] + trajs_not_nan[spawn_ix_not_nan + 1]) / 2

        if trimesh.proximity.ProximityQuery(part_mesh).signed_distance([potential_spawn]) < 0: # the spawn location is out of the part
            continue  # look for the second largest spacing (in the next loop)
        else:
            fig = plt.figure(1)
            ax = plt.gca()
            #ax = fig.add_subplot(111, projection='3d')
            x1 = [trajs_not_nan[spawn_ix_not_nan, 0], trajs_not_nan[spawn_ix_not_nan + 1, 0]]
            y1 = [trajs_not_nan[spawn_ix_not_nan, 1], trajs_not_nan[spawn_ix_not_nan + 1, 1]]
            z1 = [trajs_not_nan[spawn_ix_not_nan, 2], trajs_not_nan[spawn_ix_not_nan + 1, 2]]
            #ax.plot(x1, y1, z1, 'k')
            ax.scatter(potential_spawn[0], potential_spawn[1], potential_spawn[2], c='g', marker='o',s=4)

            trajs = np.concatenate((trajs[:spawn_ix+1], np.nan * np.zeros((1, 3, n_steps)), trajs[spawn_ix+1:]), axis=0)
            trajs[spawn_ix + 1, :, j + 1] = potential_spawn
            mass = np.concatenate((mass[:spawn_ix+1], [1], mass[spawn_ix+1:]), axis=0)
            break  # do the one spawning and stop
            
    return trajs,mass

def rem_particle(trajs,mass,j):
    """
    Remove particle at the narrowest spacing
    """

    rem_outside(trajs,mass,part_mesh,j) # first removes particles that are outside    

    for_dist = np.all(~np.isnan(trajs[:, [0, 2], j + 1]), axis=1)
    trajs_not_nan = trajs[for_dist, :, j + 1]
    trajs_not_nan_old = trajs[for_dist, :, j]
    point_dists_not_nan = compute_orth_dist(trajs_not_nan, trajs_not_nan_old)
    point_dists = np.inf * np.ones(trajs[:, :, j + 1].shape[0])
    have_diff = np.where(for_dist)[0]
    point_dists[have_diff[:-1]] = point_dists_not_nan
    kill_ix_not_nan_all = np.argsort(point_dists_not_nan)
    kill_ix_all = np.argsort(point_dists)
    
    kill_ix = kill_ix_all[0]  # kill particle at smallest distance
    
    fig = plt.figure(1)
    ax = plt.gca()
    potential_kill = trajs[kill_ix, :, j + 1]
    
    ax.scatter(potential_kill[0], potential_kill[1], potential_kill[2], c='r', marker="x",s=8,linewidth=1)  
    
    
    
    trajs[kill_ix, :, j + 1] = np.array([np.nan, np.nan, np.nan])
    mass[kill_ix] = np.nan
    
    return trajs,mass
    
def rem_outside(trajs,mass,part_mesh,j):
    """
    Checks if particles are inside or outside mesh
    Removes the particles outside
    Returns 1 if particle is outside
    """
    
    for_outside = np.all(~np.isnan(trajs[:, [0, 2], j + 1]), axis=1)
    points = trajs[for_outside, :, j + 1]
    distance = trimesh.proximity.ProximityQuery(part_mesh).signed_distance(points)
    outside_bool_not_nan = distance < 0
    
    outside_bool = np.full(trajs[:,:, j+1].shape[0], False)
    outside_bool[for_outside] = outside_bool_not_nan # vectors with True for outside particles (wrt original trajs that has NaNs too)
    
    # remove outside particles
    trajs[outside_bool, :, j + 1] = np.array([np.nan, np.nan, np.nan])
    mass[outside_bool] = np.nan
    
    return outside_bool
           

# ----------------------
# 4. Functions for plotting
# ----------------------

def plot_trajectories(trajs):
    """
    Plots trajectories in 3D using Matplotlib with interactive rotation.

    Args:
        trajs (numpy.ndarray): Array of shape (num_trajectories, num_time_steps, num_dimensions)
                              representing the trajectories to be plotted.
    """

    # Plotting
    fig = plt.figure(1)
      
    ax = plt.gca()
    #ax = fig.add_subplot(111, projection='3d')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
            
    ax.plot_trisurf(part_mesh.vertices[:, 0], part_mesh.vertices[:,1], part_mesh.vertices[:,2], triangles=part_mesh.faces, alpha = 0.1)
    
    for i in range(trajs.shape[0]):
        traj_to_plot = trajs[i, :, :]
        cmap = turbo
        c = cmap(int(i / trajs.shape[0] * 256))
        ax.plot(traj_to_plot[0, :], traj_to_plot[1, :], traj_to_plot[2, :], color=c,linewidth=0.5)

    # Add interactive rotation
    ax.view_init(elev=0, azim=90, roll=0)
    ax.set_zlim(65,85)
    ax.set_xlim(26,10)
    
    #ax.set_zlim(0,150)
    #ax.set_xlim(0,36)
    
    ax.set_aspect('equal')
    
    # Hide grid lines
    ax.grid(False)
    plt.axis('off')

    
    plt.show()
    
def plot_one_step(trajs,step):
    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax = plt.gca()
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
        
    ax.plot_trisurf(part_mesh.vertices[:, 0], part_mesh.vertices[:,1], part_mesh.vertices[:,2], triangles=part_mesh.faces, alpha = 0.2)

    for i in range(trajs.shape[0]):
        traj_to_plot = trajs[i, :, 0:step]
        cmap = turbo
        c = cmap(int(i / trajs.shape[0] * 256))
        ax.plot(traj_to_plot[0, :], traj_to_plot[1, :], traj_to_plot[2, :], color=c)

    # Add interactive rotation
    ax.view_init(elev=0, azim=90, roll=0)
    ax.set_zlim(67.5,87.5)
    ax.set_xlim(26,10)
    
    ax.set_aspect('equal')
    
        # Hide grid lines
    ax.grid(False)
    plt.axis('off')

    
    plt.savefig(cur_dir + "/" + video_frames_folder + "/frame_" + str(step).zfill(4) +'.png')
    #plt.show()
    plt.close(fig)

    

# ----------------------
# Main code
# ----------------------


# Import FEM results

mat = scipy.io.loadmat('FEM_results.mat')
part_mesh = trimesh.load_mesh('std_specimen.STL')

stress_tensors = np.array(mat['stress_tensors'])
VV = np.array(mat['VV'])
del mat

stress_tensors_swap = stress_tensors.swapaxes(0,2)

#Slow interpolation (delaunay)
#delau = scipy.spatial.Delaunay(VV) # To be computed only once! can create from existing mesh somehow
#interpolator = LinearNDInterpolator(delau, stress_tensors_swap) # Create interpolator (can be slow because of delaunay computation)

#Fast interpolation (NEAREST NEIGHBOR)
interpolator = NearestNDInterpolator(VV, stress_tensors_swap)

[PSmax, SigmaMax] = computeMaxPrincipalStress(stress_tensors)
norm_mass_factor = np.max(np.abs(SigmaMax))

# Particle steps set up

desired_distance = 0.4
n_steps = 100
step_size = 1*desired_distance # could it be changed?
K_opt = 5

# Particles beginning location
start_grid = np.transpose(np.meshgrid(np.arange(0.5, 36.5, desired_distance), [1], [61])).reshape(-1, 3)
#start_grid = np.transpose(np.meshgrid(np.arange(0.5, 36.5, desired_distance), [1], [0])).reshape(-1, 3)

useMomentum = False
inert = 0.1 # value of inertia when using momentum

trajs = np.empty((start_grid.shape[0], start_grid.shape[1], n_steps))
trajs[:] = np.NaN
trajs[:, :, 0] = start_grid


opt_duration = [] # np.zeros(n_steps)
#fig = plt.figure(1,dpi=600)
#fig.set_size_inches(3.5, 3)
fig = plt.figure(1, dpi=600)
ax = fig.add_subplot(111, projection='3d')
fig.set_size_inches(3.5, 3.5)
#fig.set_size_inches(7,7)

# repeat for each step of particle motion
   
for j in range(n_steps-1):
    
    if np.mod(j,10) == 0:
        print('Solving iteration ' + str(j))
          
    trajs_not_nan = trajs[:, :, j]
    nan_trajs_ix = np.isnan(trajs_not_nan).all(axis=1)
    trajs_not_nan = trajs_not_nan[~nan_trajs_ix, :]

    N_nodes = trajs[:, :, j].shape[0]
    mass = np.empty((N_nodes,))
    mass[:] = np.NaN

    # Compute stress at particles locations and assign particles masses
    stress = interpolateStress(trajs_not_nan, interpolator)
    [PSmax, SigmaMax] = computeMaxPrincipalStress(stress)   
    PSmax = align_psmax(PSmax)
    mass[~nan_trajs_ix] = np.abs(SigmaMax) / norm_mass_factor
    
    PSmax[np.all(~np.isnan(PSmax), axis=1), 1] = 0  # 2D version : just move in plane
    
    trajs[~nan_trajs_ix, :, j + 1] = trajs_not_nan + step_size * (PSmax)
    
    if useMomentum:
        if j > 0:
            momentum = (trajs[~nan_trajs_ix, :, j] - trajs[~nan_trajs_ix, :, j - 1])
            momentum[np.isnan(momentum)] = 0  # if it's nan (no past history) then no inertia
            momentum = momentum/np.linalg.norm(momentum) * inertia
            trajs[~nan_trajs_ix, :, j + 1] = trajs[~nan_trajs_ix, :, j + 1] + momentum
    
    #%% Optimize the position of the particles
    
    # No change in the number of particles
    trajs, cost = optimize_particles(trajs, mass, K_opt, j)  
        
    # Add one particle
    trajs_add, mass_add = add_particle(np.copy(trajs), np.copy(mass), j, part_mesh) 
    trajs_add, cost_add = optimize_particles(trajs_add, mass_add, K_opt, j)
    
    # Remove one particle
    trajs_rem, mass_rem = rem_particle(np.copy(trajs), np.copy(mass), j) #(copy.deepcopy(trajs), copy.deepcopy(mass), j)
    trajs_rem, cost_rem = optimize_particles(trajs_rem, mass_rem, K_opt, j)
    
    if cost_add < cost and cost_add < cost_rem:
        trajs = trajs_add
        mass = mass_add
        
    if cost_rem < cost and cost_rem < cost_add:
        trajs = trajs_rem
        mass = mass_rem
    # Keep the case with lowest cost
        
    # Remove particles outside of part
    rem_outside(trajs, mass, part_mesh, j)
      

#%% Plotting the results 

#fig = plt.figure(1, dpi=600)
#fig.set_size_inches(7, 7)
#ax = fig.add_subplot(111,projection='3d')

plot_trajectories(trajs)
#plt.savefig("trajectories.pdf", format="pdf", bbox_inches="tight")

#%% Compute the distances and save them

if False:
    dist = []
    for ix,traj in enumerate(trajs):
        traj = np.transpose(traj)
        other_trajs = np.delete(trajs,ix,axis=0)
        other_trajs = np.hstack(other_trajs)
        for point in traj:
            if not np.any(np.isnan(point)):
                this_dists = np.linalg.norm(np.transpose(other_trajs) - point,axis=1)
                dist.append(np.nanmin(this_dists))
    
    import pickle
    # Saving the objects:
    with open('dists.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(dist, f)


#%% Compute the stress alignment

if False:
    alignment = []
    mass = []
    for traj in trajs:
        traj = np.transpose(traj)
        for ix in range(len(traj)-1):
            if not np.any(np.isnan(traj[ix+1])) and not np.any(np.isnan(traj[ix])):
                print_dir = traj[ix+1] - traj[ix]
                print_dir = print_dir/np.linalg.norm(print_dir)
                local_stress = interpolateStress(traj[ix], interpolator)
                [PSmax, SigmaMax] = computeMaxPrincipalStress(local_stress) 
                PSmax = PSmax.squeeze()
                mass.append(np.abs(SigmaMax)/norm_mass_factor)
                alignment.append(np.abs(np.dot(print_dir,PSmax)))
       
    beta =  np.sum([a*b for a,b in zip(mass,alignment)]) /   np.sum(mass)
    print(beta)
    print(np.mean(alignment))


#%% Make a video
if False:
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    video_frames_folder = "video_frames"
    
    # Create frames storage folder
    if not os.path.exists(cur_dir + '/' + video_frames_folder):
        os.makedirs(video_frames_folder)
        
    # Delete previously stored frames in that folder    
    for filename in glob.glob(cur_dir +  '/' + video_frames_folder + '/*.png'):
        #print(filename)
        os.remove(filename)

    for i in range(n_steps):
        plot_one_step(trajs,i)
        
    # Make a video
    img_array = []
    for filename in glob.glob(cur_dir +  '/' + video_frames_folder + '/*.png'):
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        img_array.append(img)
     
    out = cv2.VideoWriter('animation.mp4',cv2.VideoWriter_fourcc(*'h264'), 5, size)
    
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()




