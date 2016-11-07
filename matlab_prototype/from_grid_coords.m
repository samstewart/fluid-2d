% converts from grid coordinates to normalized [0, 1]^2 coordinates
function [x,y] = from_grid_coords(u, v, dx)
    N = floor(1/dx);
    
    x = (u - 1)/N;
    y = (v - 1)/N;
end