% simple GUI to test advection
function test_advection()
    % more advanced: 'animating stream particles'
    
    % the main solver
    main_solver = main()
    
    % simple: just advect based on rotational velocity field
    N = 40;
    velocity_field = grid2d(N, FieldTypes.VectorField2D, @(v) [v(1) / N 0]);
    
    main_solver.dt = 1;
    total_timesteps = 100;
    
    % random grid of particles that we will advect
    particles = particle_field(N);
%     particles.values = zeros(N + 1, N + 1);
%     particles.values(N / 2, N / 2) = 1;
    
    velocity_field.plot()
    hold on;
    p = particles.plot();
    hold off;
    
    for i = 1:total_timesteps
        main_solver.advect(velocity_field, particles);
        particles.update_plot(p);
        
        drawnow limitrate;
    end
end