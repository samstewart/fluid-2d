% solves a poisson equation using jacobi method.
% We allow the user to specify two parameters that distinguish between the
% poisson equation for dissipation and the poisson equation for pressure.
% stencil: allows you to specify the coefficients of the sampled grid
% points.
%   x_{i - 1, j} stencil(1) + x_{i + 1, j} stencil(2) + x_{i, j - 1}
%   stencil(3) + x_{i, j + 1} stencil(4) + b_grid_{i, j} stencil(5)

function interior_solution = solve_poisson(b_grid, stencil, iterations) 
    % we do 35 iterations to solve the equation.
    % Method to solve a linear system via jacobi iteration
    
    % the solution will be a vector field (false)
    interior_solution = grid2d(grid.N(), false);
    
    for k = 1:iterations
        
        for i = 2:(grid.N() - 1)
            for j = 2:(grid.N() - 1)
               
               center = [i j];
               up = center + [0 1];
               down = center + [0 -1];
               right = center + [1 0];
               left = center + [-1 0];

               % now we evolve to the next step
               interior_solution.values(i, j, :) = ...
               (stencil(1) * interior_solution.gridValueAt(left) ...
               + stencil(2) * interior_solution.gridValueAt(right) ...
               + stencil(3) * interior_solution.gridValueAt(down) ...
               + stencil(4) * interior_solution.gridValueAt(up) ...
               + stencil(5) * b_grid.gridValueAt(center));

            end
        end
        
    end
end