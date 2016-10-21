% plots a 3D array representing a vector field in grid coordinates (1, 1),
% (1, 2), etc. The plotting will be in physical coordinates.
function plot_field(field)
    N = size(field, 1);
    P = size(field, 3);
    
    if P == 2
        % plot the vector field
        [x,y] = meshgrid(0:1/N:(1-1/N), 0:1/N:(1 - 1/N));
    
        % vectors transform differently than space. They don't transform
        % affinely (no translation)
        Fx = field(:, :, 2) / N;
        Fy = -field(:, :, 1) / N;
        
        % notice how we flip the coordinates to place it into physical
        % space.
        quiver(x,y, Fx, Fy);
        
    else
        % it's a scalar field so print a heatmap
        colormap('hot');
    
        imagesc(field);
    end
    
end