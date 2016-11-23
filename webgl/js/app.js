/**
The general idea:
The Fragment Shader or Pixel shader is what renders the graphics to the screen.

Tasks broken down into pieces:
Big goal:
- Write fluid solver
- Write code to solve via stepping
	- Make data types for the various fields
		- pressure texture
		- velocity texture
		- force texture
		- density texture (or particle texture)
	- Make fragment programs or pixel shaders for each operation
		- advect (velocity field)
		- advect (density or particle field)
		- dissipate
		- divergence of velocity field
		- pressure solve
		- pressure project
		- force apply
		- handle boundaries
			- add 1px boundary to crate image.

- Write code to visualize the fluid solver
	- Visualize a scalar field
		- Just write GSL shaders with color [ Done ]
		- Write GSL shader from custom texture [ Done ]

	- Visualize a vector field
		- How do they store the vecor field?
		- How do they represent the fluid?
		- How do they render multiple things on screen? Multiple textures?
			This appears to be the trick: multiple frame buffers which are double buffered. There is a similar feature in Three.JS called WebGLRenderTarget. Maybe I can construct something double buffered with this?

			You need the double buffered thing so that we don't overwrite the data as we're processing it?

			One then binds a shader to this renderer. We render a quad into each of the render targets which contains our textures

			How to do this:
			https://github.com/haxiomic/GPU-Fluid-Experiments/blob/master/Source/GPUFluid.hx
			https://github.com/haxiomic/gltoolbox/blob/master/gltoolbox/render/RenderTarget2Phase.hx

			http://webglfundamentals.org/webgl/lessons/webgl-fundamentals.html on Framebuffers

			Architecutre:
			Main code -> Fluid -> RenderTarget -> Geometry -> Shader

			We have multiple render targets and each one has a quad. The render targets can swap between two textures that they use as render targets, and the shader has access to both. The shader does the actual computation while the render targets just handle rendering. 

	- Visualize a scalar field
	- Add mouse interactivity.
Q: how do we combine multiple texture to produce the final texture?
*/
var renderer, mainRenderTarget;

var WIDTH = 400.0;
var HEIGHT = 400.0;

$(document).ready(function() {
	init();
	animate();	
});


function enableAxes() {
	// plots the axes:
	// red: X, green: Y, blue: Z
	var axisHelper = new THREE.AxisHelper( 4 );
	mainRenderTarget.scene.add(axisHelper);
}

function setupRenderer() {
	renderer = new THREE.WebGLRenderer();
	renderer.setPixelRatio( window.devicePixelRatio );
	renderer.setSize( WIDTH, HEIGHT );
	renderer.autoClear = false;

	document.body.appendChild( renderer.domElement );
}
function setupMainRenderTarget() {


	// load the texture 
	var texture = new THREE.TextureLoader().load("images/crate.jpg", function(texture) {
		texture.format = THREE.RGBAFormat;

		var shaderMaterial = new THREE.ShaderMaterial( {
					uniforms: {
						"texture": { value: texture }
					},
					vertexShader: $("#vertexshader").text(),
					fragmentShader: $("#fragmentshader").text(),
					depthWrite: false
				} );

		mainRenderTarget = new RenderTarget(shaderMaterial, WIDTH, HEIGHT);

		enableAxes();

	});
	
}


function init() {
	setupRenderer();

	setupMainRenderTarget();	
}

function animate() {
	requestAnimationFrame( animate );

	//  controls.update();

	render();		
};

// render stuff
function render() {

	if (mainRenderTarget) {

		// render the main quad with texture
		renderer.render( mainRenderTarget.scene, mainRenderTarget.camera );	
	}
	
}