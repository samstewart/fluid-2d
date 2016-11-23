

/** merges two objects together */
function mergeProperties(obj1, obj2) {

	for (var prop in obj2) {
		obj1[prop] = obj2;
	}

	return obj1;
}