mkdir build && cd build
gfortran -c ../src/env.f90 ../src/interfaces.f90 ../src/utils.f90 ../src/conjugate_gradient.f90 ../src/differential_evolution.f90 ../src/particle_swarm.f90 ../src/robust_regression.f90
ar -r optimization.a *.o && cd ../