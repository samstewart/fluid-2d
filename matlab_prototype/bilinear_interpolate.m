% simple function that interpolates the four corner points of every point into a new central vertex using bilinear interpolation. 
% field: the field that is defined only at integer grid points. This could be a vector or scalar field. Really anything with addition defined.
% coords: a 2D array of coordinates whose values we wish to compute (via interpolation) from field. Should be in grid coordinates.
% returns: if field is nxn, then the result will be (n - 1) x (n - 1) since we are computing points in the middle of every four vertices
function interpolated = bilinear_interpolate(field, coords)
	% we find the four surrounding points in the integer lattice.
	% these will be the values we use to interpolate and find the value of our point.
	% we only need the upper left grid point to construct our cell specific coordinate system.
	upper_left_grid_points = floor(coords)

	% now we find the value of the field at the four surrounding grid points for this cell.	
	upper_left_field_value = field(upper_left_grid_points)
	lower_left_field_value = field(upper_left_grid_points + [0 1])
	upper_right_grid_points = field(upper_left_grid_points + [1 0])
	lower_right_grid_points = field(upper_left_grid_points + [1 1])

	% we find our coordinates in cell space. 
	% For example, (.5, .5) would represent the middle of the cell.
	cell_coords = coords - upper_left_grid_points

	% we do bilinear interpolation (see wikipedia for exact definition).
	% The idea is that we first interpolate linearly vertically on the left and right.
	% and then we interpolate linearly horizontally.	
	
	% we first find an interpolation vertically on the left side
	left_interpolated = (cell_coords(:, :, 2) .* upper_left_grid_points + (1 - cell_coords(:, :, 2)) .* lower_left_field_value)

	% we then find an interpolate vertically on the right side.
	right_interpolated = (cell_coords(:, :, 2) .* upper_right_grid_points + (1 - cell_coords(:, :, 2)) .* lower_right_grid_points)

	interpolated = cell_coords(:, :, 1) .* left_interpolated + (1 - cell_coords(:, :, 1)) .* right_interpolated; 
end
