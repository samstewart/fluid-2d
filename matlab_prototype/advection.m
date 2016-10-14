% simple function that takes substance and transports it via a velocity field.
% currently uses an implicit method.
% velocity_field: a 2D vector array of the velocity field. [0 1/dx]^2
% substance: a 2D vector array of the substance.
% dt: change in time
% dx: size of grid cell
% returns
% advected_substance is the new array after it has been moved by the velocity field.
function advected_substance = advection(velocity_field, substance, dt, dx) 	
	N = floor(1/dx);
    advected_substance = zeros(N, N);
    
    % [TOP RIGHT BOTTOM LEFT]
    boundaries = [1 N N 1];
    
	% due to the annoyances in matlab's array notation, we cannot do this in one step but
	% must instead unwind this into a loop.
	for i = 2:(N - 1)
		for j = 2:(N - 1)
			% we move backwards in space to see where we came from.
			% we are really just following characteristics backwards
			% (looking upstream).
			grid_coords = [i j]

            % make sure we don't escape the bounds of the grid.
            % Q: does this depend on the continuity of the velocity field?
            % we have to flip the coordinates and negate the Y component to
            % fit in the array coordinates which *increase* as we go
            % downwards.
            velocity_vector_in_grid_coords = [0 -1; 1 0] * dt * N * squeeze(velocity_field(i, j, :));
            velocity_vector_in_grid_coords = velocity_vector_in_grid_coords'
			old_pos = clamp_to_range(grid_coords - velocity_vector_in_grid_coords, boundaries)
			
			% we might land "in between" grid cells so we need to interpolate.	
            bilinear_interpolate(substance, old_pos)
			advected_substance(i,j) = bilinear_interpolate(substance, old_pos);
        end
	end	
end
