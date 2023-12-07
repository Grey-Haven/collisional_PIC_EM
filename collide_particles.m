% /*uses MCC to see if the particle collided, and if so, calls CEX*/
% nn = computeNeutralDensity(particles(:,1), particles(:,2), thruster_radius);
g = sqrt(electrons(:,3).^2 + electrons(:,4).^2);

% k1 = -0.8821;
% k2 = 15.1262;
% a = k1*log(g)+k2;
% sigma = a.^2*1e-20; % collisional cross section

% mean_free_path = 1/(n_bar * hydrogen_collisional_cross_section);

sigma = pi*hydrogen_collisional_cross_section^2;

P = 1 - exp(-n_bar.*sigma.*g.*dt); % collision probability

electrons_collided_idx = find(rand(size(electrons,1),1) < P);

u_old = electrons(electrons_collided_idx,4);
v_old = electrons(electrons_collided_idx,5);

velocity_magnitude_old = sqrt(u_old.^2 + v_old.^2);

[u_new,v_new] = compute_new_velocity_burger(electrons,electrons_collided_idx);
% [u_new,v_new] = compute_new_velocity_brieda(neut_vth,electrons_collided_idx);

% alpha = rand(length(u_old),1)*2 - 1;
% u_new = alpha .* velocity_magnitude_old;
% v_new = pm(round(rand(length(u_old),1)+1)) .* sqrt(velocity_magnitude_old.^2 - u_new.^2);

assert(all(abs(sqrt(u_new.^2 + v_new.^2) - velocity_magnitude_old) < 1e-10));

% u_new = neut_vth*fmaxw(length(size(particles,1)));
% v_new = neut_vth*fmaxw(length(size(particles,1)));

electrons(electrons_collided_idx,4) = u_new;
electrons(electrons_collided_idx,5) = v_new;

% electrons(particles_collided_idx,5) = electrons(particles_collided_idx,5) + 1;

function [u_new, v_new] = compute_new_velocity_burger(particles,particles_collided_idx)

    u_old = particles(particles_collided_idx,4);
    v_old = particles(particles_collided_idx,5);

    velocity_magnitude_old = sqrt(u_old.^2 + v_old.^2);

    pm = [-1,1]';
    alpha = rand(length(u_old),1)*2 - 1;
    u_new = alpha .* velocity_magnitude_old;
    v_new = pm(round(rand(length(u_old),1)+1)) .* sqrt(velocity_magnitude_old.^2 - u_new.^2);
end

function [u_new, v_new] = compute_new_velocity_brieda(neut_vth,particles_collided_idx)
    u_new = neut_vth*fmaxw(length(particles_collided_idx));
    v_new = neut_vth*fmaxw(length(particles_collided_idx));
end

function maxwell = fmaxw(N) 
    maxwell = 2*(rand(N,1) + rand(N,1) + rand(N,1) - 1.5);
end