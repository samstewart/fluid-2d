/**
@param scale_constants an object with the various scale data (width/height/cell size). See fluid.js for full definition.
@param particleData an empty texture. We will fill it with the particle positions and velocities.
*/
function Particles(scale_constants, shaders, particleData) {

	this.particleGeometry 	= new THREE.BufferGeometry();

	// the size (in world units) of the spacing between the particles
	// we add +1 to properly center the points
	
	var particleUVs 		= new Float32Array( scale_constants.particles.total_particles * 2 );
	var positions 			= new Float32Array( scale_constants.particles.total_particles * 3 );

	// arrange the particles in model space and in texel coordintaes.
	var position = new THREE.Vector3();
	var particleUV = new THREE.Vector2();

	// fill out the particle UV coordinates
	var index = 0;

	for ( var i = 0; i < scale_constants.particles.hori_grid_points; i++ ) {
		for ( var j = 0; j < scale_constants.particles.vert_grid_points; j++ ) {

			particleUV.x = i / ( scale_constants.particles.hori_grid_points - 1 );
			particleUV.y = j / ( scale_constants.particles.vert_grid_points - 1 );

			particleUVs.toArray( particleUV, 2 * index );

			index++;
		}
	}

	// bind these attributes to the shader material

	// Note: one *MUST* add a blank position attribute to correctly initalize the mesh.
	// This behavior is undesirable since we might be setting the positions with some other method.
	this.particleGeometry.addAttribute( 'position', 	new THREE.BufferAttribute( positions, 3 ));
	this.particleGeometry.addAttribute( 'particleUV', 	new THREE.BufferAttribute( particleUVs, 2 ));

	// the particles should be white
	this.mesh = new THREE.Points( this.particleGeometry, 
		new THREE.ShaderMaterial(
		{
			uniforms: { 
				color: { value: new THREE.Color( 0xffffff ) },
				particleVariable: { value: particleData }
				size: scale_constants.particles.particle_size

			},
			fragmentShader: shaders.fragmentShader,
			vertexShader: shaders.vertexShader,
			blending: THREE.AdditiveBlending,
			transparent: false,
			depthTest: false
		})
	);
}

