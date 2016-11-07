% Helper function to convert to grid coordinates from [0, 1]^2
function [u, v] = to_grid_coords(x, y, dx)
    N = floor(1/dx);
    
    % T(x,y) = ((1 - y)(N - 1) + 1, (N - 1)x + 1)
    
    u = x * N + 1;
    v = y * N + 1;
end