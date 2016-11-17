% computes the divergence of a 2D vector field
function div = divergence(field)
   % we assume the number of rows and columns is the same.
   
   % build sparse differentiation matrix.
   % the centered finite difference stencil
   dx = 1.0 / N;
   
   
   % convert the x and y coordinate fields to vectors.
   x_stencil_field = to_stencil_vector(squeeze(field(:, :, 1)));
   y_stencil_field = to_stencil_vector(squeeze(field(:, :, 2)));
   
   
   
   % now 
   % Q: when acting on a vector, does matlab use
   N = size(field, 1);
   
   div = zeros(N, N);
end