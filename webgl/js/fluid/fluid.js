function Fluid(scale_constants, shaders, renderer) {

	this.scale_constants = scale_constants;
	this.jacobiIterations = 35;


	// a tool to manage the computations that we are doing on the shader.
	// The scale parameter controls the number of grid points in the fluid we are trying to resolve.
	this.gpuComputer = new GPUComputationRenderer(
								scale_constants.fluid.hori_grid_points, 
								scale_constants.fluid.vert_grid_points, 
								renderer);

	var pressureTexture   = this.gpuComputer.createTexture();
	var divergenceTexture = this.gpuComputer.createTexture();

	// the particle grid is actually smaller than the other textures.
	// We have fewer than one particle per point in the velocity field.
	var particleTexture	  = this.gpuComputer.createTexture(scale_constants.particles.grid_size, scale_constants.particles.grid_size);
	var velocityTexture	  = this.gpuComputer.createTexture();


	// setup the particle visualization. This will also initialize the particle texture.
	this.particles = new Particles( scale_constants, 
									shaders.particles, 
									particleTexture);

	// bind variables to textures
	// Note that we seed the particleTexture with the initial positions in the Particles constructor.
	// TODO: this separation is a bit unclear. Maybe we should do the texture initialization in here?
	this.pressureVariable 				= this.gpuComputer.addVariable('pressureVariable', shaders.pressureShader, pressureTexture);
	this.divergenceVariable 			= this.gpuComputer.addVariable('divergenceVariable', shaders.divergenceShader, divergenceTexture);
	this.particleVariable 				= this.gpuComputer.addVariable('particleVariable', shaders.particleStepShader, particleTexture);
	this.velocityVariable 				= this.gpuComputer.addVariable('velocityVariable', shaders.velocityShader, velocityTexture);

	// expose variables to the shader
	this.gpuComputer.setVariableDependencies(this.pressureVariable, [ this.pressureVariable, this.divergenceVariable ]);

	this.gpuComputer.setVariableDependencies(this.divergenceVariable, [ this.divergenceVariable ]);

	this.gpuComputer.setVariableDependencies(this.velocityVariable, [ this.velocityVariable, this.pressureVariable ]);

	this.gpuComputer.setVariableDependencies(this.particleVariable, [ this.particleVariable, this.velocityVariable ]);

	// add the extra parameters for PRESSURE shader
	this.pressureVariable.material.uniforms.dx = { value: scale_constants.fluid.dx };

	// add the extra parameters for the PARTICLE stepping shader
	this.particleVariable.material.uniforms.dt = { value: scale_constants.dt };
	this.particleVariable.material.uniforms.dragCoefficient = { value: scale_constants.particles.drag_coeff };

	// init the renderer
	var error = this.gpuComputer.init();

	if (error !== null) {
		console.error( error);
	}

	
	// center the particles properly in the screen.
	this.particles.mesh.position.z = .0001; // need to move it away from the camera a bit?
	this.particles.mesh.position.x -= scale_constants.window.width / 2.0;
	this.particles.mesh.position.y -= scale_constants.window.height / 2.0;

}

// steps the simulation by updating all of the shaders. We need a renderer to do this.
Fluid.prototype.step = function(dt) {

	// update the particle positions according to the velocity field
	this.gpuComputer.computeVariable(this.particleVariable);

	// update the actual particle mesh
	this.particles.mesh.material.uniforms.particleVariable.value = this.gpuComputer.getCurrentRenderTarget( this.particleVariable ).texture;

	// var pixelValue = new Float32Array(4 * this.scale_constants * this.particleGridSize);

	// renderer.readRenderTargetPixels(this.gpuComputer.getCurrentRenderTarget( this.particleVariable ), 0, 0, this.particleGridSize, this.particleGridSize, pixelValue);	

	// then update the actual particlel field
	// this.particles.mesh.material.uniforms.particleData.value = this.gpuComputer.getCurrentRenderTarget( this.particleVariable ).texture;

	// solve the pressure equation
	// this.solvePressure();
};

Fluid.prototype.solvePressure = function() {

	// setup the information for the pressure texture
	for (var i = 0; i < this.jacobiIterations; i++) {
	
		this.gpuComputer.computeVariable(this.pressureVariable);

		// TODO: sample a pixel value for debugging
		// var sampledPixel = new Float32Array(4);

		// this.gpuComputer.renderer.readRenderTargetPixels(this.gpuComputer.getCurrentRenderTarget( this.pressureVariable ), 256, 256, 1,1, sampledPixel);
		// console.log(sampledPixel);
	}

	// Note: the final computed pressure field will be in the read buffer
}
