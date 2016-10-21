% converts from grid coordinates to normalized [0, 1]^2 coordinates
function [x,y] = from_grid_coords(u, v, dx)
    N = floor(1/dx);
    
    x = (v - 1)/(N - 1);
    y = (1 - u)/(N - 1) + 1;
    
end