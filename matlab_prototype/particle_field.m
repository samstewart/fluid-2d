% TODO: should really extend the object grid2d and not just wrap it.
classdef particle_field < grid2d
    
    
    properties (Dependent)
        % array of coords of the particles
        particles
    end
    
    methods
        % returns N x 2 coordinates representing the particle locations in
        % the grid.
        function coords = get.particles(obj)
            
            % we make a list of nonzero coordinates
            coords = [];
            
            for i = 1:(obj.N + 1)
                for j = 1:(obj.N + 1)
                    if obj.values(i,j) > 0
                        coords = [coords; i j];
                    end
                end
            end
        end
        
        function obj = particle_field(resolution)
            % we store the particle field as a grid. 
            % 1: particle is present
            % 0: particle is not present
            obj = obj@grid2d(resolution, FieldTypes.ScalarField);
            
            % create random binary matrix
            obj.values = randi(2, obj.N + 1, obj.N + 1) - ones(obj.N + 1, obj.N + 1);
        end
        
        function update_plot(obj, p)
            coords = obj.particles();
            
            if size(coords, 1) > 0
                set(p, 'XData', coords(:, 1), 'YData', coords(:, 2));
            else
                set(p, 'XData', [], 'YData', []);
            end
        end
        
        function p = plot(obj)
            coords = obj.particles();
            p = plot(coords(:, 1), coords(:, 2), 'o', 'MarkerFaceColor', 'red');
        end
    end
end