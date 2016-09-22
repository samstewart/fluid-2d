% simple function that takes substance and transports it via a velocity field.
% currently uses an implicit method.
% velocity_field: a 2D vector array of the velocity field. [0 1/dx]^2
% substance: a 2D vector array of the substance.
% dt: change in time
% dx: size of grid cell
% 
function advected_substance = advection(velocity_field, substance, dt, dx) 	
	N = floor(1/dx);
	grid_coords = meshgrid(1:N, 1:N);

	% we move backwards in space to see where we came from.
	% we are really just following characteristics backwards
	% (looking upstream).
	old_pos = grid_coords - dt * N * velocity_field;
	
	advected_substance = quantity(old_position);
end
