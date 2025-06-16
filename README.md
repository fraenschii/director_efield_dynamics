# director_efield_dynamics
Matlab code to simulate dynamics of a model of liquid crystal director coupled to an electric field.

A finite element method to solve a coupled hyperbolic elliptic system modeling the dynamics of a liquid crystal director field subjected to an electric field. 

Run experimentj_fem.m where j=1,...,5 to run the examples in the paper. 

The file wme_fe.m is the main file, the others are subfunctions callled in it. Itcontains the actual numerical method, takes an an input a struct with the necessary parameters and initial conditions and outputs approximations of the variables d,w and \phi (mphi) at the end time of the simulation, as well as a variable 'Em' containing the reduced energy, 'damping' containing an approximation of the damping term. df, wf, phif and Ef are arrays containing approximations of the director field, the angular momentum, the electric potential and the electric field at different time steps. They can be used in combination with the m-files movie2d.m and movie2dquiver.m to create movies of the simulations. The variables tt and tf contain the time steps (tt for all time steps and tf for the ones contained in  df, wf, etc.).

movie2d.m and movie2dquiver.m can be used to create a movie with the output variables phif (movie2d.m) and df, wf, Ef (movie2dquiver.m) and time steps vector tf.

wmefield_hybrid.m is an older version of the numerical method that uses finite differences for d and w instead of finite elements.
