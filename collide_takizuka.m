


index_offset = 1;

a_x = x(1);
a_y = y(1);

lc_x = index_offset+(particles(:,1) - a_x)./dx;
lc_y = index_offset+(particles(:,2) - a_y)./dy;
is = floor(lc_x);
js = floor(lc_y);

n_L = n_bar_ions;
m_ab = (m_XE * m_XE)/(m_XE + m_XE);

m_a = m_XE;
m_b = m_XE;

coulomb_log = 10;

for i_idx = 1:Nx
    for j_idx = 1:Ny
        colliding_indices = find(is == i_idx & js == j_idx);
        particles_colliding = particles(colliding_indices);
        % shuffled_particles_colliding = particles_colliding(randperm(size(particles_colliding,1),:));

        shuffle_indices = randperm(size(particles_colliding,1));

        for k = 1:2:length(shuffle_indices)-2
            idx1 = shuffle_indices(k);
            idx2 = shuffle_indices(k);
            v1x = particles(k,3);
            v1y = particles(k,4);
            v2x = particles(k+1,3);
            v2y = particles(k+1,4);

            ux_boosted = v2x - v1x;
            uy_boosted = v2y - v1y;

            % The magnitude, another way of looking at this though is 
            % that it is rotated into a reference frame that is going only
            % in the y direction
            u_boosted = sqrt(ux_boosted^2 + uy_boosted^2);

            variance = (q_ion^2 * q_ion^2 * n_L * coulomb_log) / (8*pi*eps_0 * m_ab^2 * u_boosted^3) * dt;

            std_dev = sqrt(variance);

            R = std_dev*randn();

            assert(abs(R) <= pi/2);

            theta = 2 * atan(R);

            % ux_next = u_boosted*sin(theta);
            % uy_next = u_boosted*cos(theta);

            % dux = ux_next - ux_boosted;
            % duy = uy_next - uy_boosted;

            dux = u_boosted*sin(theta);
            duy = -u_boosted*(1 - cos(theta));

            ux_next = ux_boosted + dux;
            uy_next = uy_boosted + duy;

            v1x_next = v1x + m_ab/m_a*dux;
            v1y_next = v1y + m_ab/m_a*duy;

            v2x_next = v2x - m_ab/m_b*dux;
            v2y_next = v2y - m_ab/m_b*duy;

            v1 = norm([v1x,v1y]);
            v2 = norm([v2x,v2y]);
            v1_next = norm([v1x_next,v1y_next]);
            v2_next = norm([v2x_next,v2y_next]);

            v_prev = norm([v1x,v1y,v2x,v2y]);
            v_next = norm([v1x_next,v1y_next,v2x_next,v2y_next]);

            assert(abs((v_prev-v_next)/v_prev) < 1e-12);
            
            % assert(abs((v1-v1_next)/v1) < 1e-10);
            % assert(abs((v2-v2_next)/v2) < 1e-10);
            % assert(sqrt(v1x_next^2 + v1y_next^2 + v2x_next^2 + v2y_next) == sqrt(v1x^2 + v1y^2 + v2x^2 + v2y^2))

            particles(k,3) = v1x_next;
            particles(k,4) = v1y_next;
            particles(k+1,3) = v2x_next;
            particles(k+1,4) = v2y_next;

        end

    end
end