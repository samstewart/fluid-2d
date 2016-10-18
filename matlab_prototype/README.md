A simple fluid simulator based on the wonderful article [1].

1. Harris, M. "Fast Fluid Dynamics Simulation on the GPU," *GPU Gems*. Available at [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html).


Where I left off:
I am making a convenience function to generate a force field in user coordinates [0, 1]^2. This will allow me to test the advection code for different fields.
For example, those with swirl, etc.f