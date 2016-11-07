% plots a 3D array representing a vector field in grid coordinates (1, 1),
% (1, 2), etc. The plotting will be in grid coordinates so the field
% should be in grid coordinates.
function plot_field(field)
    N = size(field, 1) - 1; % we always have one less dimension.
    P = size(field, 3);
    
    if P == 2
        % plot the vector field
        [x,y] = meshgrid(1:1:(N + 1), 1:1:(N + 1));
   
        % we have to adjust the coordinates since now according to quiver()
        % we have (i, j) representing (y, x)
        Fx = transpose(field(:, :, 1));
        Fy = transpose(field(:, :, 2))

        % note the parameter of zero to prevent scaling.
        quiver(x,y, Fx, Fy);        
    else
        % it's a scalar field so print a heatmap
        colormap('hot');
    
        % we transpose the field so that we have the column and row
        % orderings corresponding to our grid construction.
        imagesc(flipud(transpose(field)));
    end
end