#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#*******************************
# Reset memory
#*******************************
from IPython import get_ipython
get_ipython().run_line_magic('reset', '-sf') 

import numpy as np
import matplotlib.pyplot as plt
import math
from ase.visualize import view
from ase.io import read
from numpy.linalg import inv
from ase.build import make_supercell
from ase import Atoms
from numba import njit
import pandas as pd



def read_cell_file(datafile):
    """
    Function that reads the cell file.
    
    Parameters:
        datafile: String
            Path of the file.
            
    returns: np.ndarray
        Contains all the values of the cell file.
    """
    
    with open(datafile, "r") as file:
        cell = file.read()
        cell = cell.split("\n")
        cell = np.array([c.split() for c in cell[1:-1]]).astype(float)
        return cell



@njit
def Sf_histogram(grid, q_val, super_pos, particles):
    """
    Function that creates a structure factor histogram for a single q-value,
    where the X-axis represents the distance and the Y-axis represents 
    the structure factor values distributed over distance intervals.

    Parameters:
        grid: np.ndarray
            The X-axis values of the histogram, representing atomic pair distances.
        q_val: int
            A single q-value appearing in the structure factor.
        super_pos: np.ndarray
            The positions of all particles in the supercell, represented as vectors.
        particles: int
            The number of particles in the supercell.

    Returns: np.ndarray
        The structure factor histogram for a single q-value, with values assigned 
        to their corresponding distance intervals.
    """
    
    grid_dr = grid[1]-grid[0]
    
    #Number of particles in the supercell
    N = super_pos.shape[0]
    
    q_D = super_pos*q_val
    histo = np.zeros_like(grid)
    
    
    for i in range(0, particles):
        for j in range(0,N):
            if not i == j:
                norm = 0.0
                for ic in range(3):
                    # Distance between two particles with pythagorean theorem
                    norm += (q_D[i,ic]-q_D[j,ic])**2
                    
                dis = np.sqrt(norm)
                
                #Calculates one index value for the structure value
                increment = (1/particles)*np.sin(dis)/(dis)
                if int((dis/q_val)/grid_dr) > len(histo):
                    raise Exception(f"Gridin suurin arvo liian pieni: kahden hiukkasen välinen etäisyys iteraatiossa on {int(dis/q_val)}")
                  
                histo[int((dis/q_val)/grid_dr)] += increment
                
    return histo
            




def structure_factor(traj, cell, supercell_size, q_values, snap):
    """
    Function that constructs a supercell from the given trajectory for a single snapshot
    and computes the structure factor histogram for all q-values from the supercell.

    Parameters:
        traj: np.ndarray
            Trajectory data.
        cell: np.ndarray
            Cell data.
        supercell_size: int
            The size of the supercell, e.g., a value of 3 creates a 3x3x3 supercell.
            The number must be odd to ensure that the original cell is at the center of the supercell.
        q_values: list
            A list of q-values on which the structure factor depends.
        snap: int
            An integer indicating the snapshot number of the trajectory.

    Returns: np.ndarray, np.ndarray
        The first np.ndarray is the grid variable from the histogram used in the weighting function,
        and the second np.ndarray contains the structure factor histograms for each q-value.
    """
    
    
    if supercell_size % 2 == 0:
        raise ValueError("The size of the supercell is not odd")
        

    traj_frame = traj[snap]
    cell_frame = cell[snap]
    
    #cell vectors
    a = cell_frame[2:5]
    b = cell_frame[5:8]
    c = cell_frame[8:11]


    s_max_val = int((supercell_size+1)/2)
    s_min_val = int((1-supercell_size)/2)

    particles = len(traj_frame)
    supercell = traj_frame.positions.copy()
    
    #Adds cells to the supercell
    for i in range(s_min_val,s_max_val):
        for j in range(s_min_val,s_max_val):
            for k in range(s_min_val,s_max_val):
                if i == 0 and j == 0 and k == 0:
                    continue
                else:
                    supercell = np.vstack((supercell,traj_frame.positions + a*i+b*j+c*k))
    

    cell_diagonal = supercell_size*(12+12+12)
    grid = np.arange(0, cell_diagonal/2, 0.1)
    
    
    sf_sum = []
    for qi in q_values:
        print(f"q-val: {qi}")
        histo = Sf_histogram(grid, qi, supercell, particles)
        sf_sum.append(histo)
    
    sf_sum = np.array(sf_sum)
    return sf_sum,grid
    



def weight_function(grid: list):
    """
    Function that creates a weight function used to modify the structure factor values.

    Parameters:
        grid: np.ndarray
        This uses the same grid variable as the one returned by the structure_factor function.

    Returns: np.ndarray
        An array of the weighting function's values corresponding to the distances in the grid variable.
    """
    
    weight = []
    for r in grid:

        f = math.exp(-0.5*(r+0.05)**2/100**2)
        weight.append(f)

    weight = np.array(weight)
    return weight


"""
In this script, the structure factor for atoms is computed for multiple snapshots,
and the average of these is displayed in a plot.
Here, the structure factor is derived from a modified Debye scattering equation.
Instead of including every atomic pair present in the supercell, only the pairs
where at least one atom belongs to the original cell are considered.
The original cell is expanded into a supercell so that the diffraction peaks are
clearly distinguishable in the plot.

A weight function is also applied to the structure factor. For each q-value, a
histogram is computed that shows the distribution of the structure factor over
different distances for that particular q-value. Each histogram is then multiplied
by the weight function, modifying the significance of the atomic pair distances;
the contribution from distant pairs is nearly nullified in the structure factor,
while that of the nearby pairs is enhanced. This ensures that atoms at the edges
of the supercell do not affect the result.
"""





traj = read("npt_coords.xyz", index=":") 
cell = read_cell_file("/npt_cell.txt")

expt_e = pd.read_csv("/epsilon_mittaus.csv",header = None)
expt_d = pd.read_csv("/delta_data.csv",header = None)


#Changes the experimental data that is compared to the simulation data
pressure = 10


supercell_size = 11
q_values = np.linspace(2.0,3.4,100)


#Calculates the structure factor for each snapshot
aver_sf = []
for snap in range(0,len(cell)):
    print(f"snapshot: {snap}")
    sf_sum,grid = structure_factor(traj,cell, supercell_size, q_values, snap)
    aver_sf.append(sf_sum)
aver_sf = np.array(aver_sf)



#Average is taken from each histogram and this is multiplied by the weight function
sf_plot = 1 + np.sum(weight_function(grid) * np.average(aver_sf,axis=0),axis = 1)
sf_plot_max = np.max(sf_plot)


# Plotting the graph of the average structure factor
plt.figure()
plt.plot(q_values,sf_plot,label= "Sim.")
if pressure == 10:
    plt.plot(expt_d[0][0:61],expt_d[1][0:61]*(sf_plot_max/np.max(expt_d[1])),color="black",label="Expt.")  #delta
else:
    plt.plot(expt_e[0][0:61],expt_e[1][0:61]*(sf_plot_max/np.max(expt_e[1])),color="black",label="Expt.")  #epsilon
plt.legend()
plt.grid()
plt.xlabel(r"$q$ [Å$^{-1}]$")
plt.ylabel("$S(q)$ [arb. units]")
plt.title(f"{pressure} GPa")

