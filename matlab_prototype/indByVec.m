% allows you to index a multidimensional array by a vector
function elem = indByVec(array, vec)
	elem = array(vec2ind(size(array), vec));
end
