A simple fluid simulator based on the wonderful article [1].

1. Harris, M. "Fast Fluid Dynamics Simulation on the GPU," *GPU Gems*. Available at [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html).


Where I left off:
There is a problem with plotting vector fields of the form [0 y]. They should face straight up, but they aren't.

Next step, compute the pressure field by using jacobi iteration to solve a Poisson equation

\nabla^2 P = \nabla w (velocity field with divergence included)