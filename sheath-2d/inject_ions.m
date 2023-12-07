function particles = inject_ions(N_particles, tag, charge, mass, weight, ...
                                 x_min, x_max, y_min, y_max)
    particle_locations = zeros(N_particles,3);
    particle_velocities = zeros(N_particles,3);
    particle_tags = zeros(N_particles,1);
    particle_charges = zeros(N_particles,1);
    particle_masses = zeros(N_particles,1);
    particle_weights = zeros(N_particles,1);
%     del_x = 1/(sqrt(N_particles)-1);
%     del_y = 1/(sqrt(N_particles)-1);
    avg_v_mag = 1/mass;
    
    x_samples = x_min + (x_max-x_min) .* rand(N_particles,1);
    y_samples = y_min + (y_max-y_min) .* rand(N_particles,1);

    try
        for idx = 1:N_particles
                particle_locations(idx,:) = [x_samples(idx), y_samples(idx), 0];
                particle_velocities(idx,:) = avg_v_mag*[randn(), randn(), 0];
                particle_tags(idx) = tag;
                particle_charges(idx) = charge;
                particle_masses(idx) = mass;
                particle_weights(idx) = weight;
        end
        particles = [particle_locations,particle_velocities,particle_tags,particle_masses,particle_weights,particle_charges];
    catch
        disp("foo");
    end
end