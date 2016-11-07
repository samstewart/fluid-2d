% creates a 2D field (in grid coordinates) based on the two helper
% functions that you pass in.
% F: at F(x, y) you should the vector at that point. x, y will be in [0,
% 1]^2.
% dx: the resolution for the field
function field = create_field(F, dx)
    N = floor(1/dx);
    
    velocity_field = zeros(N, N, 2);
    
    for i = 1:(N + 1)
        for j = 1:(N + 1)
            % put the coordinate in [0, 1]^2.
            [x,y] = from_grid_coords(i, j, dx);
            
            v = F([x y]);
            
            % vectors transform under pushforward.
            u = v(1) * N;
            v = v(2) * N;
            
            % now put into coordinates.
            velocity_field(i, j, :) = [u v];
        end
    end
    
    field = velocity_field;
end