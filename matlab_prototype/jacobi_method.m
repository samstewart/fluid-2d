% Method to solve a linear system via jacobi iteration
% A: matrix in Ax = b
% b: column vector in Ax = b
% N: number of iterations
% returns: column vector solution after N iterations

function sol = jacobi_method(A, b, N)
	diagonal = diag(diag(A)); % strip out the diagonal
	diag_deleted = A - diagonal; % delete the diagonal
	
	sol = zeros(size(b, 1), 1); % initial guess of zero
    
	for i = 1:N
        % computing the matrix inverse
		sol = diagonal \ (b - diag_deleted * sol)
    end
    
end
