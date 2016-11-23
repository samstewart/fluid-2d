% represents 2D grid of region [1, N]x[1,N] (we call these grid
% coordinates)
classdef grid2d < handle
    
    
    properties
        % 2D array of values (i, j) |-> (x_i, y_i)
        values 
        % number of grid points
        N
        % the type of field (scalar or vector field)
        field_type
    end
    
    properties (Dependent)
        % stencil for computing pressure poisson equation
        pressure_stencil
        % stencil for computing dissipation equation for velocity
        dissipation_stencil
        % dirichlet boundary conditions for velocity
        velocity_boundary_stencil
        % neumann boundary conditions for pressure
        pressure_boundary_stencil
        
        % spatial resolution
        dx
        % size of the grid (N + 1) x (N + 1)
        grid_size
        % the divergence of this vector field
        divergence
    end
    
    methods (Access = private)
        
        
        % simple command that converts a vector to a multidimensional index given the size of the array
        % you should ignore any singleton dimensions (don't provide an element of the vector)
        % TODO turn this into a Gist
        % size_of_matrix: vector of dimensions of the matrix (can include singleton dimensions)
        % vec_index: a vector whose nth component corresponds to the nth index. 
        function indx = vec2ind(obj, size_of_matrix, vec_index)
            % We convert this to a linear index. Not that matlab takes the full n dimensional integer
                % lattice and converts it to a one dimensional integer lattice by reading off successive columns.	
            % the formula (really bijection between countable sets) between Z^n_1 x Z^n_2 x Z^n_3 --- Z^n_k is given by 
            % (a_0, a_1, ..., a_n) -> a_0 + (b_0)(a_1 - 1) + (b_0 b_1) (a_2 - 1) + (b_0 b_1 b_2) (a_3 - 1) + ....
            % article idea: multidimensional array indexing and showing that cross products are countable.

            % find the multiplying factors (how many columns to skip)

            % filter out singleton dimensions
            size_no_singleton = size_of_matrix(size_of_matrix > 1);

            % note that we don't need the final product since this would place us off by one.
            multipliers = [1 cumprod(size_no_singleton(1:end - 1))];

            % subtract 1 from all but the first index
            vec_index(2:end) = vec_index(2:end) - 1;

            % then we compute the formula to convert from n dimensions down to one dimension.
            indx = sum(multipliers .* vec_index);
        end
        
        % clamps the given coordinate to the range passed.
        % coord: the given pair of coordinates
        % 
        function clamped = clamp_to_boundaries(obj, coord)
            % range: [TOP RIGHT BOTTOM LEFT]
            range = [(obj.N + 1) (obj.N + 1) 1 1];
            clamped(1) = max(range(4), min(range(2), coord(1)));
            clamped(2) = max(range(3), min(range(1), coord(2)));
        end
        
    end
    
    methods
        function set.N(obj, N_new) 
            obj.N  = N_new;
        end
        
        function dx = get.dx(obj)
            dx = 1.0 / obj.N;
        end
        
        function N = get.N(obj)
            N = obj.N;
        end
        
        function g_size = get.grid_size(obj)
            s      = size(obj.values);
            g_size = s(1 : end - 1);
        end    
        % stencil for computing pressure poisson equation
        function s = get.pressure_stencil(obj)
            % you can derive this from the symbolic definition of Jacobi
            % iteration
            s = 1/4 * [1 1 1 1 -obj.dx()^2];
        end
        
        % stencil for computing dissipation equation for velocity
        function s = get.dissipation_stencil(obj)
            s = zeros(1, 4);
        end
        % dirichlet boundary conditions for velocity
        function s = get.velocity_boundary_stencil(obj)
            s = -obj.pressure_boundary_stencil();
        end
        
        function s = get.pressure_boundary_stencil(obj)
            s = [1 1 1 1];
        end
        
        % compute the divergence of this grid using centered finite
        % differences. Does nothing at the boundaries; these must be
        % handled separately.
        function div = get.divergence(obj)
            % fill a grid with zeros
            div = grid2d(obj.N, false);
            
            % handle all non-boundary values
            for i = 2:obj.N()
                for j = 2:obj.N()
                    % update the x-coordinate
                    div.values(i, j, 1) = 1 / (2 * obj.dx) * (obj.values(i + 1, j) - obj.values(i - 1, j));
                    
                    % update the y-coordinate
                    div.values(i, j, 2) = 1 / (2 * obj.dx) * (obj.values(i, j + 1) - obj.values(i, j - 1));
                end
            end
        end
        % plots a 3D array representing a vector field in grid coordinates (1, 1),
        % (1, 2), etc. The plotting will be in grid coordinates so the field
        % should be in grid coordinates.
        function p = plot(obj)
            if obj.field_type == FieldTypes.ScalarField
                % it's a scalar field so print a heatmap
                colormap('hot');

                % we transpose the field so that we have the column and row
                % orderings corresponding to our grid construction.
                p = imagesc(flipud(transpose(obj.values)));
                
            else
                % plot the vector field
                [x,y] = meshgrid(1:1:(obj.N + 1), 1:1:(obj.N + 1));

                % we have to adjust the coordinates since now according to quiver()
                % we have (i, j) representing (y, x)
                Fx = transpose(obj.values(:, :, 1));
                Fy = transpose(obj.values(:, :, 2));

                % note the parameter of zero to prevent scaling.
                p = quiver(x,y, Fx, Fy);        
            end
        end
        
        % updates a given plot depending on our field type
        function update_plot(obj, p)
            if obj.field_type == FieldTypes.ScalarField
                

                % we transpose the field so that we have the column and row
                % orderings corresponding to our grid construction.
                p.CDataMapping = 'scaled';
                p.CData = 
                
            else
                % TODO: update vector field plot       
            end
        end
    end
        
    methods (Access = public)
        
        
        % constructs a grid with N sample points in each direction (will
        % actually be N + 1 since we include zero).
        % resolution: first argument specifies the grid resolution. If < 1
        % then we assume its dx. If > 1 then we assume its N.
        % type: field type from FieldTypes enum.
        % F: a function for constructing the vector field dynamically. For
        % example, F([x y]) = v.
        function obj = grid2d(varargin)
            % Set the resolution of the grid
            resolution = varargin{1};
            
            if resolution < 1
                obj.N = floor( 1 / resolution );
            else
                obj.N = resolution;
            end

            % decide what kind of field to construct
            obj.field_type = varargin{2};
            
            if obj.field_type == FieldTypes.ScalarField
                % scalar field
                obj.values = zeros(obj.N + 1, obj.N + 1, 1);
            else
                % 2D vector field
                obj.values = zeros(obj.N + 1, obj.N + 1, 2);
            end
            
            
            % if we are given an additional function for constructing the
            % vector field, we use it.
            if nargin > 2
                % we were given a function for constructing the velocity
                % field
                F = varargin{3};
                
                for i = 1:(obj.N + 1)
                    for j = 1:(obj.N + 1)
                        % put the coordinate in [0, 1]^2.
                        [x,y] = obj.from_grid_coords(i,j);

                        v = F([x y]);

                        % vectors transform under pushforward so only the scaling.
                        u = v(1) * obj.N;
                        v = v(2) * obj.N;

                        % now put into coordinates.
                        obj.values(i, j, :) = [u v];
                    end
                end
            end

        end
        
        % converts from grid coordinates to normalized [0, 1]^2 coordinates
        function [x,y] = from_grid_coords(obj, u, v)
            x = (u - 1) / obj.N;
            y = (v - 1) / obj.N;
        end

        % Helper function to convert to grid coordinates from [0, 1]^2
        function [u, v] = to_grid_coords(obj, x, y)
            % T(x,y) = ((1 - y)(N - 1) + 1, (N - 1)x + 1)
            u = x * obj.N + 1;
            v = y * obj.N + 1;
        end

        % takes a scalar field and fills out all the values from the points already given
        % using bilinear interpolation. Creates a new MxM interpolated field.
        function refine_grid(obj, M)
            
            interpolated = zeros(M, M);

            % now we go through and interpolate each entry based on the four
            % surrounding entries. 
            for i = 1:(M + 1)
                for j = 1:(M + 1)
                    % fill in the grid value using bilinear interpolation
                    interpolated(i, j) = obj.interpolateValueAt([1 + (i - 1)/M 1 + (j - 1)/M]);
                end
            end
            
            obj.values = interpolated;
        end

        % gets the grid value at the given vector coordinate. Clamps
        % coordinate to boundaries.
        function value = gridValueAt(obj, coord)
            coord = obj.clamp_to_boundaries(coord);
            
            value = squeeze(obj.values(coord(1), coord(2), :));
        end
        
        % boundary_stencil: [RIGHT TOP LEFT BOTTOM]
        function set_boundaries(obj, boundary_stencil)
            
            % handle the right boundary (little clunky, but slice notation
            % is clunky too).
            % Note: I'm not totally sure how to handle the corners, but now
            % just update each of them twice.
            for i = 1:(obj.N + 1)
                obj.values(end, i) = boundary_stencil(1) * obj.values(end - 1, i);
            end
            
            % handle the top boundary
            for i = 1:(obj.N + 1)
                obj.values(i, end) = boundary_stencil(2) * obj.values(i, end - 1);
            end
            
            % handle the left boundary
            for i = 1:(obj.N + 1)
                obj.values(1, i) = boundary_stencil(3) * obj.values(2, i);
            end
            
            % handle the bottom boundary
            for i = 1:(obj.N + 1)
                obj.values(i, 1) = boundary_stencil(4) * obj.values(i, 2);
            end

        end
        % simple function that interpolates the four corner points of every point into a new central vertex using bilinear interpolation. 
        % Interpolates the four vertices surrounding the coordinate given by
        % 'coords'.
        % field: the field that is defined only at integer grid points. This could be a vector or scalar field. Really anything with addition defined.
        % coord: a point in grid coordinates where (0, 0) is upper left and (N, N)
        % is lower right. The first coordinate is the row, and the second
        % coordinate is the column. This is essentially "array index" coordinates
        % and is the most natural for internal storage.
        % returns: will return the interpolated value. This type depends on the type of the field (which should be a two dimensional vector field).
        function interpolated = interpolateValueAt(obj, coord)
            
            % TODO: shorten code by writing everything as matrices.
            
            % we find the four surrounding points in the integer lattice.
            % these will be the values we use to interpolate and find the value of our point.
            % we only need the upper left grid point to construct our cell specific coordinate system.
            lower_left_grid_point = floor(coord);

            % now we find the value of the field at the four surrounding grid points for this cell.	
            lower_left_field_value   = obj.gridValueAt(lower_left_grid_point);
            lower_right_field_value  = obj.gridValueAt(lower_left_grid_point + [1 0]);
            upper_right_field_value  = obj.gridValueAt(lower_left_grid_point + [1 1]);
            upper_left_field_value   = obj.gridValueAt(lower_left_grid_point + [0 1]);

            % we find our coordinates in cell space. 
            % For example, (.5, .5) would represent the middle of the cell.
            cell_coords = coord - lower_left_grid_point;

            % we intentionally flip the coordinates since the row represents the y
            % coordinate and the column represents the x coordinate.
            cell_x = cell_coords(1);
            cell_y = cell_coords(2);

            % we do bilinear interpolation (see wikipedia for exact definition).
            % The idea is that we first interpolate linearly vertically on the left and right.
            % and then we interpolate linearly horizontally.	

            % we first find an interpolation vertically on the left side
            left_interpolated = (1 - cell_y) * lower_left_field_value + cell_y * upper_left_field_value;

            % we then find an interpolate vertically on the right side.
            right_interpolated = (1 - cell_y) * lower_right_field_value + cell_y * upper_right_field_value;

            % we now interpolate horizontally between these values
            interpolated  	   = (1 - cell_x) * left_interpolated + cell_x * right_interpolated;
            
        end
    end
end
