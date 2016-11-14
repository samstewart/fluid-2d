% solves the poisson equation:
% \nabla^2 p = \nabla \cdot w
% to correct the new vector field w after one step to a vector field with
% zero divergence. This method then subtracts \nabla p from w to find a
% field that is divergence free.
% One can view the pressure p as a Lagrange multiplier in
% this sense. 
% Interesting perspective on Helmholtz-Hodge decomposition: http://www2.cs.uh.edu/~chengu/Teaching/Spring2013/Lecs/Lec12.pdf
function projected_field = pressure_project(original_field, dx, dt)
    
end