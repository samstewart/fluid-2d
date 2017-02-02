/** 
Fills a given THREE.FloatType texture with colors given by the map function.
map (grid_coord, uv_coord) 
	grid_coord = {
		i: [ the horizontal index],
		j: [ the vertical index ],
		index: the linear index of the current pixel
	},
	uv_coord = {
		u: [ horizontal texel coordinate],
		v: [ vertical texel coordinate ]
	}

coordinates are such that origin is in lower left corner.

2/1/2017: this method is currently deprecated.
*/
function fillTexture(texture, map) {
	var data = texture.image.data;

	// note: we need to divide by four.
	for (var i = 0; i < data.length / 4; i++ ) {
		// ordering believed to be from:
		// https://github.com/mrdoob/three.js/blob/6c7f000734f8579da37fb39e5c2e9e5e2dfb14f8/examples/js/utils/ImageUtils.js
		var grid_coord = {
			i: i % texture.image.width, 
			j: Math.floor(i / texture.image.width),
			index: i
		}

		var uv_coord = {
			u: grid_coord.i / (texture.image.width - 1),
			v: grid_coord.j / (texture.image.height - 1),
		}

		var pixelColor = map(grid_coord, uv_coord);

		data[i * 4] 	= pixelColor[0];
		data[i * 4 + 1] = pixelColor[1];
		data[i * 4 + 2] = pixelColor[2];
		data[i * 4 + 3] = pixelColor[3];
	}

	texture.needsUpdate = true; // not completely sure if this is necessary?
}