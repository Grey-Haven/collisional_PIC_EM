function [preserved_particles,preserved_prev_particle_locations] = remove_particles(particles, prev_particle_locations, x_min, x_max, y_min, y_max)    
    particle_locations = particles(:,1:3);
    preserved_particle_idxs = find(particle_in_zone(particle_locations(:,:), x_min, x_max, y_min, y_max));
    preserved_particles = particles(preserved_particle_idxs,:);
    preserved_prev_particle_locations = prev_particle_locations(preserved_particle_idxs,:);
end