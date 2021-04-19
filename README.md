# director_efield_dynamics
Matlab code to simulate dynamics of a model of liquid crystal director coupled to an electric field.

A mixed finite difference/finite element method to solve a coupled hyperbolic elliptic system modeling the dynamics of a liquid crystal director field subjected to an electric field. 

The file wmefield_hybrid.m contains the actual numerical method, takes an an input a struct with the necessary parameters and initial conditions and outputs approximations of the variables d,w and \phi (mphi) at the end time of the simulation, as well as a variable 'Em' containing the reduced energy, 'damping' containing an approximation of the damping term, 'maxgrad' containing the evolution of \max_x\|\Grad_h d_h\| over time. df, wf, phif and Ef are arrays containing approximations of the director field, the angular momentum, the electric potential and the electric field at different time steps. They can be used in combination with the m-files movie2d.m and movie2dquiver.m to create movies of the simulations. The variables tt and tf contain the time steps (tt for all time steps and tf for the ones contained in  df, wf, etc.).

The files num_exp1.m and num_exp2.m generate structs that can be used as input data for wmefield_hybrid.m, initdata1.mat can be loaded in the workspace and also contains a struct that can be used as input data for wmefield_hybrid.m.

movie2d.m and movie2dquiver.m can be used to create a movie with the output variables phif (movie2d.m) and df, wf, Ef (movie2dquiver.m) and time steps vector tf.
