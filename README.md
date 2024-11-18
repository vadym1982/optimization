gfortran -c env.f90 utils.f90 interfaces.f90 test_functions.f90 particle_swarm.f90 differential_evolution.f90 conjugate_gradient.f90 robust_regression.f90
ar -r optimization.a *.o