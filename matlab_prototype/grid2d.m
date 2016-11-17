% represents 2D grid of region [1, N]x[1,N] (we call these grid
% coordinates)
classdef grid2d < handle
    properties
        % 2D array of values (i, j) |-> (x_i, y_i)
        values 
        
        % discretized laplacian operator (we cache the matrix)
        discretized_laplacian_operator
        
        % discretized divergence operator (we cache the matrix)
        discreteized_divergence_operator
        
        % number of grid points
        N
        % determines if scalar or vector field
        is_scalar_field
    end
    
    properties (Dependent)
        % this grid represented as a stencil
        stencil_vector
        % size of the stencil vector
        stencil_length
        
        % spatial resolution
        dx
        % size of the grid (N + 1) x (N + 1)
        grid_size
    end
    
    methods (Access = private)
        % creates a sparse matrix from a stencil for four points. 
        % The stencil should have the following ordering when converting
        % from the 2D grid to 1D vector
        %       4
        %       |
        %       |
        % 2 --- 3 ----- 5
        %       |
        %       |
        %       1
        
        function m = create_matrix_from_stencil(stencil)
            % make bands from stencil into sparse matrix

            % first band
            i1 = 1:N;
            j1 = 1:N;
            s1 = ones(1, N) * stencil(1);

            % second band
            j2 = i1 + 1;
            s2 = ones(1, N) * stencil(2);

            % third band
            j3 = i1 + 2;
            s3 = ones(1, N) * stencil(3);

            % fourth band
            j4 = i1 + 3;
            s4 = ones(1, N) * stencil(4);

            % fifth band
            j5 = i1 + 4;
            s5 = ones(1, N) * stencil(5);

            m = sparse([i1 i1 i1 i1 i1], [j1 j2 j3 j4 j5], [s1 s2 s3 s4 s5]);
        end
        
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
        
        function stencil_length = get.stencil_length(obj)
            stencil_length = 5 * (obj.N - 1)^2;
        end
        
        % code that converts the vector we use when multiplying by the stencil
        % operator back into a grid. The ordering we are "unconverting" is
        %       4
        %       |
        %       |
        % 2 --- 3 ----- 5
        %       |
        %       |
        %       1

        function set.stencil_vector(obj, vector_in)
            obj.values = zeros(obj.N, obj.N);

            % the current index of the center point in above diagram in the
            % stencil_vector.
            cur_center_point_index = 3;
            
            % we stay away from the boundaries because they don't have all the
            % neighbors.
            % now we update the grid.
            for i = 2:(obj.N - 1)
                for j = 2:(obj.N - 1)
                    % convert the grid to a vector using the ordering described
                    % above

                    % Point 1:
                    obj.values(i, j - 1) = vector_in(cur_center_point_index - 2);

                    % Point 2:
                    obj.values(i - 1, j) = vector_in(cur_center_point_index - 1);

                    % Point 3:
                    obj.values(i, j) = vector_in(cur_center_point_index);

                    % Point 4:
                    obj.values(i, j + 1) = vector_in(cur_center_point_index + 1);

                    % Point 5:
                    obj.values(i + 1, j) = vector_in(cur_center_point_index + 2);

                    % now advance the index by 6 so that we advance to the next
                    % center node
                    cur_center_point_index = cur_center_point_index + 5; 
                end
            end
        end
        
        function is_scalar_field = get.is_scalar_field(obj) 
            is_scalar_field = (size(obj.values, 3) == 1);
        end
        
        % plots a 3D array representing a vector field in grid coordinates (1, 1),
        % (1, 2), etc. The plotting will be in grid coordinates so the field
        % should be in grid coordinates.
        function plot(obj)
            if obj.is_scalar_field()
                % it's a scalar field so print a heatmap
                colormap('hot');

                % we transpose the field so that we have the column and row
                % orderings corresponding to our grid construction.
                imagesc(flipud(transpose(obj.values)));
                
            else
                % plot the vector field
                [x,y] = meshgrid(1:1:(obj.N + 1), 1:1:(obj.N + 1));

                % we have to adjust the coordinates since now according to quiver()
                % we have (i, j) representing (y, x)
                Fx = transpose(obj.values(:, :, 1));
                Fy = transpose(obj.values(:, :, 2));

                % note the parameter of zero to prevent scaling.
                quiver(x,y, Fx, Fy);        
            end
        end

        % converts a scalar grid to stencil ordering in one huge vector so that we
        % can discretize the linear operators.
        % The ordering is as follows:
        %       4
        %       |
        %       |
        % 2 --- 3 ----- 5
        %       |
        %       |
        %       1
        function stencil_vector = get.stencil_vector(obj)

            % Q: how much more inefficient is it to allow the repetition of
            % neighbors during this conversion process?

            % total_rows = 5 points per interior grid point
            stencil_vector = zeros( obj.stencil_length(), 1);

            % the current index of the center point in above diagram in the
            % stencil_vector.
            cur_center_point_index = 3;

            % we stay away from the boundaries because they don't have all the
            % neighbors.
            for i = 2:(obj.N - 1)
                for j = 2:(obj.N - 1)
                    % convert the grid to a vector using the ordering described
                    % above

                    % Point 1:
                    stencil_vector(cur_center_point_index - 2)  = field(i, j - 1);

                    % Point 2:
                    stencil_vector(cur_center_point_index - 1)  = field(i - 1, j);

                    % Point 3:
                    stencil_vector(cur_center_point_index)      = field(i, j);

                    % Point 4:
                    stencil_vector(cur_center_point_index + 1)  = field(i, j + 1);

                    % Point 5:
                    stencil_vector(cur_center_point_index + 2)  = field(i + 1, j);

                    % now advance the index by 6 so that we don't overwrite the old
                    % values
                    cur_center_point_index = cur_center_point_index + 5; 
                end
            end
        end
        
    end

    methods (Access = public)
        % constructs a grid with N sample points in each direction (will
        % actually be N + 1 since we include zero).
        % resolution: first argument specifies the grid resolution. If < 1
        % then we assume its dx. If > 1 then we assume its N.
        % is_scalar: specifies if it is a scalar field or not.
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
            is_scalar = varargin{2};
            if is_scalar
                obj.values = zeros(obj.N + 1, obj.N + 1, 1);
            else
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
            
            % pre-compute the discretized differential operators
            % now build the differentiation operator from the centered finite
            % difference stencil
            stencil = 1/(2*obj.dx)* [ 0 -1 0 1 0];
            
            % compute the divergence operator
            obj.discreteized_divergence_operator = create_matrix_from_stencil(stencil, obj.stencil_length());
            
            stencil = 1/(2 * dx) * [1 1 -4 1 1];
            
            % compute the laplacian operator
            obj.discreteized_laplacian_operator  = create_matrix_from_stencil(stencil, obj.stencil_length());
            
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
            
            % we index into the values array by converting the object into
            % a vector.
            value = squeeze(obj.values(obj.vec2ind(size(obj.values), coord)));
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
