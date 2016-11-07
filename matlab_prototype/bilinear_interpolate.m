% simple function that interpolates the four corner points of every point into a new central vertex using bilinear interpolation. 
% Interpolates the four vertices surrounding the coordinate given by
% 'coords'.
% field: the field that is defined only at integer grid points. This could be a vector or scalar field. Really anything with addition defined.
% coord: a point in grid coordinates where (0, 0) is upper left and (N, N)
% is lower right. The first coordinate is the row, and the second
% coordinate is the column. This is essentially "array index" coordinates
% and is the most natural for internal storage.
% returns: will return the interpolated value. This type depends on the type of the field (which should be a two dimensional vector field).
% TIME GOES BY TOO QUICKLY.
function interpolated = bilinear_interpolate(field, coord)
	% we find the four surrounding points in the integer lattice.
	% these will be the values we use to interpolate and find the value of our point.
	% we only need the upper left grid point to construct our cell specific coordinate system.
	lower_left_grid_point = floor(coord);
    
    % we need to make sure we don't move off the boundary
    N = size(field, 1) - 1;
    
    % [TOP RIGHT BOTTOM LEFT]
    boundary = [(N + 1) (N + 1) 1 1];
    
	% now we find the value of the field at the four surrounding grid points for this cell.	
	lower_left_field_value   = indByVec(field, lower_left_grid_point);
	lower_right_field_value  = indByVec(field, clamp_to_range(lower_left_grid_point + [1 0], boundary));
	upper_right_field_value  = indByVec(field, clamp_to_range(lower_left_grid_point + [1 1], boundary));
	upper_left_field_value   = indByVec(field, clamp_to_range(lower_left_grid_point + [0 1], boundary));

	% we find our coordinates in cell space. 
	% For example, (.5, .5) would represent the middle of the cell.
	cell_coords = coord - lower_left_grid_point
    lower_left_grid_point
    lower_right_grid_point =  clamp_to_range(lower_left_grid_point + [1 0], boundary)
    lower_left_field_value
    lower_right_field_value
    
    % we intentionally flip the coordinates since the row represents the y
    % coordinate and the column represents the x coordinate.
	cell_x = cell_coords(1);
    cell_y = cell_coords(2);
	
	% we do bilinear interpolation (see wikipedia for exact definition).
	% The idea is that we first interpolate linearly vertically on the left and right.
	% and then we interpolate linearly horizontally.	
	
	% we first find an interpolation vertically on the left side
	left_interpolated = (1 - cell_y) * lower_left_field_value + cell_y * upper_left_field_value

	% we then find an interpolate vertically on the right side.
	right_interpolated = (1 - cell_y) * lower_right_field_value + cell_y * upper_right_field_value
	
	% we now interpolate horizontally between these values
	interpolated  	   = (1 - cell_x) * left_interpolated + cell_x * right_interpolated
end
