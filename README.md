# GMX_Advanced_Tools
This repository stores tools to extract customized data from molecular dynamics trajectories enabling parallelization via CUDA and MPI.
The tools are developed as the forks to the original GMX HYDORDER tool. 
The patches enable compilation togather and common interface of command line options with Gromacs v.2025.
Legacy C-style coding of the original tool is preserved for smooth integration.

The tools included in this repository:
HYDORDER2
Calculates distribution of the Orientational Tetrahonal Order (OTO) parameter. This is similar to the original HYDRODER tool, but the output is histogram of the OTO calculated within the specified region (all box or within (or outside) a shell around "Protein" group, the shell is defined by lower and upper radius). Furthermore a more convinient rescaled form of OTO parameter is added, providing the result normalized in [0:1] interval. See the tool individual description for more details on usage.

BARRELFIT
This tool fits the approximate cyllindrical surface to the could of points representing the heavy atoms of protein's backbone. The current implementation assumes the protein to have a barrel-like closed surface and fits it with two cylinderssmoothely connected with a sphere.  
