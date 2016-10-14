A simple fluid simulator based on the wonderful article [1].

1. Harris, M. "Fast Fluid Dynamics Simulation on the GPU," *GPU Gems*. Available at [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html).


Where I left off:
I am trying to write the function *interpolate_field* to create a fully interpolated field from just four corner values. I'm having some scaling problems.

I also need to make the code consistent so that we use the [row column] coordinates for internal representation. I don't want to use too many conversion functions.

I am having trouble with the advection code when we try for a horizontal velocity field. The vertical advection used to work, but now I changed coordinates.