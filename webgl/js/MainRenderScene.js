// simple class to manage the scene we use for rendering
function MainRenderScene(width, height) {
	this.renderer = renderer;

	this.width = width;
	this.height = height;

	////////////////////////////////////////////////////
	// setup the main scene
	////////////////////////////////////////////////////
	this.scene = new THREE.Scene();

	// create the lights one lighting the front
	var light = new THREE.DirectionalLight(0xffffff);
	light.position.set(0, 0, 1);
	this.scene.add(light);

	// create the camera
	this.camera = new THREE.OrthographicCamera(- this.width / 2.0, this.width / 2.0, this.height / 2.0, -this.height / 2.0, -1000, 1000);

	// place the rendering quad in the scene

	this.renderQuad = this.createRenderQuad();
	//this.scene.add(this.renderQuad);
}


MainRenderScene.prototype.getPassThroughVertexShader = function() {
	return	"varying vec2 vUv;\n"
			"void main()	{\n" +
			"\n" +
			"\nvUv = uv;\n"+
			"	gl_Position = vec4( position, 1.0 );\n" +
			"\n" +
			"}\n";

}

MainRenderScene.prototype.getPassThroughFragmentShader = function( w, h ) {

	return  "\nvarying vec2 vUv;\n"
			"\nuniform sampler2D outputTexture;\n\n" +
			"void main() {\n" +
			"\n" +
			"	gl_FragColor = texture2D(texture, vUv);\n" +
			"\n" +
			"}\n";

}

/** constructs a unit quad quad with the given shader (which might include a texture) */
MainRenderScene.prototype.createRenderQuad = function() {

	// visualize the solution to the pressure equation.
	var shaderMaterial = new THREE.ShaderMaterial( {
				uniforms: {
					"outputTexture": { value: null }
				},
				vertexShader:   this.getPassThroughVertexShader(),
				fragmentShader: this.getPassThroughFragmentShader(),
				depthWrite: false
			} );

	var plane = new THREE.PlaneGeometry( this.width, this.height );

	// we use a shader to render the geometry
	var quad = new THREE.Mesh( plane, shaderMaterial );
	quad.matrixAutoUpdate = false;

	return quad;

}

