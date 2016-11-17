% simple function that takes substance and transports it via a velocity field.
% currently uses an implicit method.
% velocity_field: a 2D vector array of the velocity field. Should be specified in grid coords.
% substance: a 2D vector array of the substance.
% dt: change in time
% dx: size of grid cell
% returns
% advected_substance is the new array after it has been moved by the velocity field.
function advected_substance = advection(velocity_grid, substance_grid, dt) 	
    % we need to allow for boundaries.
    advected_substance = zeros(substance_grid.grid_size());
 
	% due to the annoyances in matlab's array notation, we cannot do this in one step but
	% must instead unwind this into a loop.
    % we update only the interior cells.
	for i = 2:substance_grid.N
        for j = 2:substance_grid.N
			% we move backwards in space to see where we came from.
			% we are really just following characteristics backwards
			% (looking upstream).
            
			grid_coords = [i j]
            
            % we use the stable scheme by looking upstream. 
			old_pos = grid_coords - dt * squeeze(velocity_field.values(i, j, :))';
			
			% we might land "in between" grid cells so we need to interpolate.	
			advected_substance(i,j) = substance_grid.interpolateValueAt(old_pos);
        end
	end	
end