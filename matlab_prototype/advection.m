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


	% due to the annoyances in matlab's array notation, we cannot do this in one step but
	% must instead unwind this into a loop.
	for i = 2:(N - 1)
		for j = 2:(N - 1)
			% we move backwards in space to see where we came from.
			% we are really just following characteristics backwards
			% (looking upstream).
			grid_coords = [i j];

			old_pos = grid_coords - dt * N * velocity_field(i, j);
			
			% we might land "in between" grid cells so we need to interpolate.	
			advected_substance(i,j) = bilinear_interpolate(substance, old_pos);
		end
	end	
end
