% takes a scalar field and fills out all the values from the points already given
% using bilinear interpolation. Creates a new NxN interpolated field.
function interpolated = interpolate_field(field, N)
    interpolated = zeros(N, N);
    
    % now we go through and interpolate each entry based on the four
    % surrounding entries. 
    for i = 1:(N + 1)
        for j = 1:(N + 1)
            interpolated(i, j) = bilinear_interpolate(field, [1 + (i - 1)/N 1 + (j - 1)/N]);
        end
    end
    
end