% main code that solves fluid equation
classdef main
    properties
        % the size of the grid
        N = 10
        % the size of the timestep
        dt = .1;
    end
    
    methods
        
        % steps the solution by applying the following:
        % 1. advect the fluid
        % 2. dissipate the velocity field due to viscosity
        % 3. apply the forces
        % 4. project the fluid back to the kernel of the divergence
        % operator.
        function new_u = step_solution(obj, u)
            u = grid2d(obj.N(), FieldTypes.VectorField2D);
            
            % you first advect the velocity field
            dt = 1; % TODO: need to recover this from the ODE solver.
            u = obj.advect(u, u);
            
            % then we let the viscosity dissipate the vector field
            u = obj.dissipate(u);
            
            % then we add the forces
            u = obj.add_forces(u);
            
            % then we project back onto the subspace (kernel of the
            % divergence)
            u = obj.pressure_project(u, dt);
        end
        
        % solves the poisson equation:
        % \nabla^2 p = \nabla \cdot w
        % to correct the new vector field w after one step to a vector field with
        % zero divergence. This method then subtracts \nabla p from w to find a
        % field that is divergence free.
        % One can view the pressure p as a Lagrange multiplier in
        % this sense. 
        % Interesting perspective on Helmholtz-Hodge decomposition: http://www2.cs.uh.edu/~chengu/Teaching/Spring2013/Lecs/Lec12.pdf
        function projected_field = pressure_project(obj, velocity)
            % find pressure
            % TODO: do we normalize by dt or does it matter?
            pressure = solve_poisson(velocity.divergence(), velocity.pressure_stencil(), 35); 

            % project back to the kernel of the divergence operator
            projected_field = velocity.values - pressure;
        end
        
        % simple function that takes substance and transports it via a velocity field.
        % currently uses an implicit method.
        % velocity_field: a 2D vector array of the velocity field. Should be specified in grid coords.
        % substance: a 2D vector array of the substance.
        % dt: change in time
        % dx: size of grid cell
        % returns
        % advected_substance is the new array after it has been moved by the velocity field.
        function advect(obj, velocity_grid, substance_grid) 	
            % we need to allow for boundaries.
            advected_values = zeros(substance_grid.grid_size());

            % due to the annoyances in matlab's array notation, we cannot do this in one step but
            % must instead unwind this into a loop.
            % we update only the interior cells.
            for i = 2:substance_grid.N
                for j = 2:substance_grid.N
                    % we move backwards in space to see where we came from.
                    % we are really just following characteristics backwards
                    % (looking upstream).

                    grid_coords = [i j];

                    % we use the stable scheme by looking upstream. 
                    old_pos = grid_coords - obj.dt * velocity_grid.gridValueAt(grid_coords)';

                    % we might land "in between" grid cells so we need to interpolate.	
                    advected_values(i,j) = substance_grid.interpolateValueAt(old_pos);
                end
            end	
            
            % update the grid values
            substance_grid.values = advected_values;
        end
    end
end