% solves a poisson equation 
% \nabla^2 u = b(x)
% on a rectangular domain, but does not handle the boundary conditios
function interior_solution = solve_poisson(field, b, dx)
    % build sparse matrix representing the Laplacian
    
    % convert the field the special ordering necessary for inverting the
    % matrix.
    stencil_vector = to_stencil_vector(field);
    
    % Create stencil using this ordering:
    %       4
    %       |
    %       |
    % 2 --- 3 ----- 5
    %       |
    %       |
    %       1
    stencil = 1/(2 * dx) * [1 1 -4 1 1];
    laplacian_matrix = create_matrix_from_stencil(stencil, size(stencil_vector, 1));
    
    % we do 35 iterations to solve the equation.
    interior_solution = jacobi_method(laplacian_matrix, b, 35);
end