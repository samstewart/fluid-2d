% plots a 3D array representing a vector field in grid coordinates (1, 1),
% (1, 2), etc.
function plot_field(field)
    N = size(field, 1);
    M = size(field, 2);
    P = size(field, 3);
    
    if P == 2
        % plot the vector field
        [x,y] = meshgrid(1:1:N, 1:1:M);
    
        u = field(:, :, 1);
        v = field(:, :, 2);
    
        quiver(x,y,u,v);
    else
        % it's a scalar field so print a heatmap
        colormap('hot');
    
        imagesc(field);
    end
    
end