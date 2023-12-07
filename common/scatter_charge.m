try
    index_offset = 1;
    rho = zeros(Nx+1,Ny+1);
    for p = 1:length(particles)
        particle = particles{p};
        lc_x = index_offset+(particle.location(1)+offset_x)/dx;
        lc_y = index_offset+(particle.location(2)+offset_y)/dy;
        
        % Jx_stagger_x Jx_stagger_y
        
        i = floor(lc_x);
        j = floor(lc_y);

        q = particle.charge;
        weight = weightMap(particle.tag);
        wx = q*weight;
        
        x_node = x(i) + Jx_stagger_x;
        y_node = y(j) + Jx_stagger_y;
        fx = particle.location(1) - x_node;
        fy = particle.location(2) - y_node;

        rho(i,j) = rho(i,j) + (1-fx)*(1-fy)*wx;
        rho(i,j+1) = rho(i,j+1) + (1-fx)*fy*wx;
        rho(i+1,j) = rho(i+1,j) + fx*(1-fy)*wx;
        rho(i+1,j+1) = rho(i+1,j+1) + fx*fy*wx;
    end
catch exception
    throw(exception);
end