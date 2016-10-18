% creates a 2D field (in grid coordinates) based on the two helper
% functions that you pass in.
% F: at F(x, y) you should the vector at that point. x, y will be in [0,
% 1]^2.
% dx: the resolution for the field
function field = create_field(F, dx)
    N = floor(1/dx);
    
    field = zeros(N, N);
    
    for i = 1:N
        for j = 1:N
            % put the coordinate in [0, 1]^2.
            user_space_coord = [1 - (i / N) (j / N)];
            [x,y] = F(user_space_coord);
            field(i, j) = [N* ;
            
        end
    end
end