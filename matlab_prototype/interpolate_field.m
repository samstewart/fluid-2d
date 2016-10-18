% takes a scalar field and fills out all the values from the points already given
% using bilinear interpolation. Creates a new M*N interpolated field.
function interpolated = interpolate_field(field, M, N)
    interpolated = zeros(N, M);
    
    % now we go through and interpolate each entry based on the four
    % surrounding entries. Of course we must skip the boundary values. 
    for i = 2:(N - 1)
        for j = 2:(M - 1)
            interpolated(i, j) = bilinear_interpolate(field, [1 + i/N 1 + j/M]);
        end
    end
    
end