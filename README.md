# director_efield_dynamics
Matlab code to simulate dynamics of a model of liquid crystal director coupled to an electric field.

A mixed finite difference/finite element method to solve a coupled hyperbolic elliptic system modeling the dynamics of a liquid crystal director field subjected to an electric field. 

The file wmefield_hybrid.m contains the actual numerical method, takes an an input a struct with the necessary parameters and initial conditions and outputs approximations of the variables d,w and $$\phi$$ at later times. 
