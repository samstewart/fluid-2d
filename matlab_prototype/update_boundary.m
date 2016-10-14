% updates a given field by pulling values from one layer away from the
% boundary.
% field: the actual field we wish to update.
% coefficients: a vector of four numbers that allows us to update both the
% pressure and the velocity fields (we'll need different coefficients for
% each)
%   For velocity field [-1 -1 -1 -1]
%   For the pressure field [1 1 1 1] (because we are using finite
%   difference)
%   1: the coefficient for the left boundary
%   2: the coefficient for the right boundary
%   3: the coefficient for the top boundary
%   4: the coefficient for the bottom boundary
% returns the field after updating the boundary
function updated_field = update_boundary(field, coefficients)
	% the actualy boundary lies between two cells. The no slip condition
	% gives that, for example on the left,
    %   \frac{u_0,j + u_1,j}{2} = 0
    % which gives the condition
    %   u_0,j = -u_1, j
    
    % handle the left boundary
    field(:, 1, :) = coefficients(1) * field(:, 2, :);
    
    % handle the right boundary
    field(:, end, :) = coefficients(2) * field(:, end - 1, :);
    
    % handle the top boundary
    field(1, :, :) = coefficients(3) * field(2, :, :);
    
    % handle the bottom boundary
    field(end, :, :) = coefficients(4) * field(end - 1, :, :);
    
    % TODO: does it matter what we set the four corners to?
    updated_field = field;
end
