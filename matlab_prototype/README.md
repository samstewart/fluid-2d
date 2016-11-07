A simple fluid simulator based on the wonderful article [1].

1. Harris, M. "Fast Fluid Dynamics Simulation on the GPU," *GPU Gems*. Available at [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html).


Where I left off:

TODO:
1. Convert bilinear interpolation and advection code to use new coordinate system.
2. Write Poisson equation solver
3. Convert whole representation to a staggered grid.

Next step, compute the pressure field by using jacobi iteration to solve a Poisson equation

\nabla^2 P = \nabla w (velocity field with divergence included)

