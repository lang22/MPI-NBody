# MPI-NBody
MPI version of the Barnes Hut algorithm, coded for C with GLFW for visualization

The code should be compiled right away with the makefile. Pay special attention to dependencies (GLFW has a LOT of them).
Once compiled, only requires 3 parameters to execute:

mpiexec -np 2 ./NBODY.x <no_of_Particles> <no_of_Steps> <write_whatever_to_activate_liveview>

so for example

mpiexec -np 2 ./NBODY.x 2000 10 i-like-pineaples
will execute the liveview with 2000 particles (in this case, the number of steps is ignored, it runs non-stop)
(And the third argument doesn't matter, it is just to get the argument count to 4)
