% simple command that converts a vector to a multidimensional index given the size of the array
% you should ignore any singleton dimensions (don't provide an element of the vector)
% TODO turn this into a Gist
% size_of_matrix: vector of dimensions of the matrix (can include singleton dimensions)
% vec_index: a vector whose nth component corresponds to the nth index. 
function indx = vec2ind(size_of_matrix, vec_index)
	% We convert this to a linear index. Not that matlab takes the full n dimensional integer
        % lattice and converts it to a one dimensional integer lattice by reading off successive columns.	
	% the formula (really bijection between countable sets) between Z^n_1 x Z^n_2 x Z^n_3 --- Z^n_k is given by 
	% (a_0, a_1, ..., a_n) -> a_0 + (b_0)(a_1 - 1) + (b_0 b_1) (a_2 - 1) + (b_0 b_1 b_2) (a_3 - 1) + ....
	% article idea: multidimensional array indexing and showing that cross products are countable.

	% find the multiplying factors (how many columns to skip)

	% filter out singleton dimensions
	size_no_singleton = size_of_matrix(size_of_matrix > 1);

	% note that we don't need the final product since this would place us off by one.
	multipliers = [1 cumprod(size_no_singleton(1:end - 1))];
	
	% subtract 1 from all but the first index
	vec_index(2:end) = vec_index(2:end) - 1;
	
	% then we compute the formula to convert from n dimensions down to one dimension.
	indx = sum(multipliers .* vec_index);
end
