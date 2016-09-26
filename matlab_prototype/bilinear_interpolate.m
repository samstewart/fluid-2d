% simple function that interpolates the four corner points of every point into a new central vertex using bilinear interpolation. 
% field: the field that is defined only at integer grid points. This could be a vector or scalar field. Really anything with addition defined.
% coords: a point in grid coordinates that we wish to interpolate. 
% returns: will return the interpolated value. This type depends on the type of the field (which should be a two dimensional vector field.
function interpolated = bilinear_interpolate(field, coord)
	% we find the four surrounding points in the integer lattice.
	% these will be the values we use to interpolate and find the value of our point.
	% we only need the upper left grid point to construct our cell specific coordinate system.
	upper_left_grid_point = floor(coord)

	% now we find the value of the field at the four surrounding grid points for this cell.	
	upper_left_field_value  = indByVec(field, upper_left_grid_point); 
	lower_left_field_value  = indByVec(field, upper_left_grid_point + [0 1]); 
	upper_right_field_value = indByVec(field, upper_left_grid_point + [1 0]);
	lower_right_field_value  = indByVec(field, upper_left_grid_point + [1 1]);

	% we find our coordinates in cell space. 
	% For example, (.5, .5) would represent the middle of the cell.
	cell_coords = coord - upper_left_grid_point
	cell_x = cell_coords(1)
	cell_y = cell_coords(2)
	
	% we do bilinear interpolation (see wikipedia for exact definition).
	% The idea is that we first interpolate linearly vertically on the left and right.
	% and then we interpolate linearly horizontally.	
	
	% we first find an interpolation vertically on the left side
	left_interpolated = cell_y * upper_left_field_value + (1 - cell_y) * lower_left_field_value

	% we then find an interpolate vertically on the right side.
	right_interpolated = cell_y * upper_right_field_value + (1 - cell_y) * lower_right_field_value
	
	% we now interpolate horizontally between these values
	interpolated  	   = cell_x * left_interpolated + (1 - cell_x) * right_interpolated;
end
