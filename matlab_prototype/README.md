A simple fluid simulator based on the wonderful article [1].

1. Harris, M. "Fast Fluid Dynamics Simulation on the GPU," *GPU Gems*. Available at [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html).


Where I left off:

I am making a convenience function to generate a force field in user coordinates [0, 1]^2. This will allow me to test the advection code for different fields.
For example, those with swirl, etc.

As I just learned in a graphics talk, I need to make a connection between grid coordinates and the unit square. 

This will be an affine transformation that I use whenever plotting or inputing data.

There is a problem with plotting vector fields of the form [0 y]. They should face straight up, but they aren't.

Next step, compute the pressure field by using jacobi iteration to solve a Poisson equation

\nabla^2 P = \nabla w (velocity field with divergence included)

