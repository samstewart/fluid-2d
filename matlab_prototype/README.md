A simple fluid simulator based on the wonderful article [1].

1. Harris, M. "Fast Fluid Dynamics Simulation on the GPU," *GPU Gems*. Available at [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html).




TODO:
1. Finish code for converting the grid to a stencil vector. Convert divergence calculation and poisson solver to sparse matrix calculation. 
2. Write Poisson equation solver
3. Convert whole representation to a staggered grid.

Next step, compute the pressure field by using jacobi iteration to solve a Poisson equation

\nabla^2 P = \nabla w (velocity field with divergence included)

