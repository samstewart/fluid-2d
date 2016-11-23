% Solves a poisson equation using Jacobi method.
% We allow the user to specify two parameters that distinguish between the
% poisson equation for dissipation and the poisson equation for pressure.
% stencil: allows you to specify the coefficients of the sampled grid
% points.
%   x_{i - 1, j} stencil(1) + x_{i + 1, j} stencil(2) + x_{i, j - 1}
%   stencil(3) + x_{i, j + 1} stencil(4) + b_grid_{i, j} stencil(5)
% stencil for pressure equation:
% [1/4 1/4 1/4 1/4 -b_grid.dx()^2 / 4]
% 
% stencil for dissipation equation:
%   not sure how he's choosing viscosity. I think he's choosing \delta t v
%   = 1? But an n appears out of nowhere.
%
% boundary_stencil: the coefficients for constructing the boundary.
% 
% [ RIGHT TOP LEFT BOTTOM ]
% 
% Boundary stencil for velocity field. We have the Neumann boundary
% condition u = 0 on the boundary. Since the boundary line lies between two
% grid points, we say
% (u_{1, N - 1} + u_{1, N})/2 = 0
% which produces a constraint on the velocity field.
% -[1 1 1 1]
% Boundary stencil for pressure field:
% [1 1 1 1]

% Test cases:
% 1: you can test this against a known solution on (-1, 1)x(-1, 1) with
% Dirichlet conditions u = 0 on the boundary. 
% Matlab code for testing it:
% u = (1 - x^2)/2 - 
%    16/\[Pi]^3 Sum[
%      Sin[k \[Pi] ( 1 + x ) /2 ]/( k^3 Sinh[k \[Pi]] )* ( 
%        Sinh[k \[Pi] (1 + y)/2] + Sinh[k \[Pi] (1 - y )/2]), {k, 1, 
%       10, 2}];
% v = u /. {x -> 2 s - 1, y -> 2 r - 1};
% DensityPlot[v, {r, 0, 1}, {s, 0, 1}]

function interior_solution = solve_poisson(b_grid, stencil, iterations, boundary_stencil) 
    % we do 35 iterations to solve the equation.
    % Method to solve a linear system via jacobi iteration
    
    % the solution will be a vector field (false)
    interior_solution = grid2d(b_grid.N(), b_grid.field_type);
    
     % TODO: strange bug where we get the negative of the solution we
     % should be seeing.
    for k = 1:iterations
        for i = 2:b_grid.N()
            for j = 2:b_grid.N()
               
               center = [i j];
               up = center + [0 1];
               down = center + [0 -1];
               right = center + [1 0];
               left = center + [-1 0];

               % now we evolve to the next step by averaging over the
               % adjacent grid cells.
               interior_solution.values(i, j, :) = ...
               (stencil(1) * interior_solution.gridValueAt(left) ...
               + stencil(2) * interior_solution.gridValueAt(right) ...
               + stencil(3) * interior_solution.gridValueAt(down) ...
               + stencil(4) * interior_solution.gridValueAt(up) ...
               + stencil(5) * b_grid.gridValueAt(center));

            end
        end
        
        % fix the boundaries
        interior_solution.set_boundaries(boundary_stencil);
    end
    
    
    
end