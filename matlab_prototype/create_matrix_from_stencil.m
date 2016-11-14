% creates a sparse matrix from a stencil for four points. The stencil should have
% the following ordering:
%       4
%       |
%       |
% 2 --- 3 ----- 5
%       |
%       |
%       1
function m = create_matrix_from_stencil(stencil, N)
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