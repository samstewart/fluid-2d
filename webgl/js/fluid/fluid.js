function Fluid(N, renderer) {

	// we need a context in which to render
	this.renderer = renderer;
	this.N = N;
	this.jacobiIterations = 30;

	// the background buffer for rendering the pressure texture
	var pressureBuffer;

	this.prototype.solvePressure = solverPressure() {
		// setup the information for the pressure texture
		
		for (var i = 0; i < this.jacobiIterations; i++) {
			
		}
	}
}
