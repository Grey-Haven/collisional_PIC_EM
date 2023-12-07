function particles = inject_particles(N_particles, ...
                                      tag, charge, mass, weight, ...
                                      x_min, x_max, y_min, y_max)
    
%     particles = Particle.empty(N_particles, 0);
    particle_locations = zeros(N_particles,3);
    particle_velocities = zeros(N_particles,3);
    particle_tags = zeros(N_particles,1);
    particle_charges = zeros(N_particles,1);
    particle_masses = zeros(N_particles,1);
    particle_weights = zeros(N_particles,1);
    
    x_samples = x_min + (x_max-x_min) .* rand(N_particles,1);
    y_samples = y_min + (y_max-y_min) .* rand(N_particles,1);
    
    avg_v_mag = 1;
    
    for i = 1:N_particles
%         particles{i} = Particle([x_samples(i), y_samples(i), 0], ...
%                                 [v_injection, 0, 0], ... 
%                                 tag, charge, mass);
        particle_locations(i,:) = [x_samples(i), y_samples(i), 0];
        particle_velocities(i,:) = avg_v_mag*[randn(), randn(), 0];
        particle_tags(i) = tag;
        particle_charges(i) = charge;
        particle_masses(i) = mass;
        particle_weights(i) = weight;
    end
    particles = [particle_locations,particle_velocities,particle_tags,particle_masses,particle_weights,particle_charges];
end