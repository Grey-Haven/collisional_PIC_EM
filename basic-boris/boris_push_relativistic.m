Ex_padded = zeros(size(Ex,1)+2,size(Ex,2)+2,size(Ex,3)+2);
Ex_padded(2:end-1,2:end-1,2:end-1) = Ex;

Ey_padded = zeros(size(Ey,1)+2,size(Ey,2)+2,size(Ey,3)+2);
Ey_padded(2:end-1,2:end-1,2:end-1) = Ey;

Ez_padded = zeros(size(Ez,1)+2,size(Ez,2)+2,size(Ez,3)+2);
Ez_padded(2:end-1,2:end-1,2:end-1) = Ez;

Hx_padded = zeros(size(Hx,1)+2,size(Hx,2)+2,size(Hx,3)+2);
Hx_padded(2:end-1,2:end-1,2:end-1) = Hx;

Hy_padded = zeros(size(Hy,1)+2,size(Hy,2)+2,size(Hy,3)+2);
Hy_padded(2:end-1,2:end-1,2:end-1) = Hy;

Hz_padded = zeros(size(Hz,1)+2,size(Hz,2)+2,size(Hz,3)+2);
Hz_padded(2:end-1,2:end-1,2:end-1) = Hz;

% Note that we're padding and shifting x and y, just padding z
x_padded = (x(1)-dx:dx:x(end)-dx)+dx/2;
y_padded = (y(1)-dy:dy:y(end)-dy)+dy/2;
z_padded = z(1)-dz:dz:z(end)+dz;

particle_locations = particles(:,1:3);
particle_velocities = particles(:,4:6);
% particle_tags = particles(3);
particle_masses = particles(:,8);
% particle_weights = particles(:,9);
particle_charges = particles(:,10);

% particles(1:end-injection_rate)
% particles(end-injection_rate:injection_rate)

% particles_in_domain = particles(1:end-injection_rate,:);
% particles_out_domain = particles(end-injection_rate+1:end,:);



particles_in_domain_idxs = find(particle_in_zone(particles(:,:), a_x, b_x, a_y, b_y));
particles_out_domain_idxs = find(~particle_in_zone(particles(:,:), a_x, b_x, a_y, b_y));

particles_in_domain = particles(particles_in_domain_idxs,:);
particles_out_domain = particles(particles_out_domain_idxs,:);

if (size(particles_in_domain,1) ~= 0)

    xe = particles_in_domain(:,1:3);
    ve = particles_in_domain(:,4:6);
    me = particles_in_domain(:,8);
    qe = particles_in_domain(:,10);

    E = evalE(Ex_padded,Ey_padded,Ez_padded,xe,...
              x,y,z,x_padded,y_padded,z_padded,...
              dx,dy,dz,...
              Jx_stagger_x,Jx_stagger_y,0,...
              Jy_stagger_x,Jy_stagger_y,0,...
              0,0,0,... % no Ez field, no staggering
              a_x,a_y,a_z,b_x,b_y,b_z);
    B = evalB(Hx_padded,Hy_padded,Hz_padded,xe,...
              x,y,z,x_padded,y_padded,z_padded,...
              dx,dy,dz,...
              0,0,0,... % no Hx field
              0,0,0,... % no Hy field
              Hz_stagger_x,Hz_stagger_y,0,...
              a_x,a_y,a_z,b_x,b_y,b_z);


    % gamma = sqrt(1/(1 - (v/c)^2)) = sqrt(1 + (u/c)^2)

    % Convert from v to u (u = gamma*v)
    gamma = computeGammaV(ve,kappa);
    ue = gamma.*ve;

    % Compute the next velocity in terms of u
    ue_next = updateVelocity(xe,ue,qe,me,E,B,kappa,dt);

    % Convert from u to v (v = u/gamma)
    gamma = computeGammaU(ue_next,kappa);
    ve_next = ue_next./gamma;

    particles_in_domain(:,4:6) = ve_next;
    newLocation = pushParticle(xe,ve_next,dt);
    particles_in_domain(:,1:3) = newLocation;
end
if (size(particles_out_domain,1) ~= 0)

    xe = particles_out_domain(:,1:3);
    ve = particles_out_domain(:,4:6);
    me = particles_out_domain(:,8);
    qe = particles_out_domain(:,10);

%     gamma = computeGammaV(ve,kappa);
% 
%     ue = gamma*ve;
    newLocation = pushParticle(xe,ve,dt);
    particles_out_domain(:,1:3) = newLocation;
end

particles = [particles_in_domain; particles_out_domain];


function E = evalE(Ex_padded,Ey_padded,Ez_padded,loc,x,y,z,...
                   x_padded,y_padded,z_padded,dx,dy,dz,...
                   Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                   Ey_stagger_x,Ey_stagger_y,Ey_stagger_z,...
                   Ez_stagger_x,Ez_stagger_y,Ez_stagger_z,...
                   a_x,a_y,a_z,b_x,b_y,b_z)
    
    
    E1 = gather_field(Ex_padded,loc,x_padded,y,z_padded,dx,dy,dz,...
                      Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                      0,0,1,... % padding offset on x,y,z, respectively
                      a_x,a_y,a_z,b_x,b_y,b_z); % Not adding padding to y, this is due to boundaries, not padding
    E2 = gather_field(Ey_padded,loc,x,y_padded,z_padded,dx,dy,dz,...
                      Ey_stagger_x,Ey_stagger_y,Ey_stagger_z,...
                      0,0,1,... 
                      a_x,a_y,a_z,b_x,b_y,b_z);
%     E3 = gather_field(Ez_padded,loc,x,y,z_padded,dx,dy,dz,...
%                       Ez_stagger_x,Ez_stagger_y,Ez_stagger_z,...
%                       0,0,0,...
%                       a_x,a_y,a_z,b_x,b_y,b_z);
    E3 = zeros(size(E1));
    E = [E1,E2,E3];
end

function B = evalB(Hx_padded,Hy_padded,Hz_padded,loc,x,y,z,...
                   x_padded,y_padded,z_padded,dx,dy,dz,...
                   Hx_stagger_x,Hx_stagger_y,Hx_stagger_z,...
                   Hy_stagger_x,Hy_stagger_y,Hy_stagger_z,...
                   Hz_stagger_x,Hz_stagger_y,Hz_stagger_z,...
                   a_x,a_y,a_z,b_x,b_y,b_z)
    
%     B1 = gather_field(Hx_padded,loc,x,y,z,dx,dy,dz,...
%                       Hx_stagger_x,Hx_stagger_y,Hx_stagger_z,...
%                       0,0,0,...
%                       a_x,a_y,a_z,b_x,b_y,b_z);
%     B2 = gather_field(Hy_padded,loc,x,y,z,dx,dy,dz,...
%                       Hy_stagger_x,Hy_stagger_y,Hy_stagger_z,...
%                       0,0,0,...
%                       a_x,a_y,a_z,b_x,b_y,b_z);
    B3 = gather_field(Hz_padded,loc,x_padded,y_padded,z_padded,dx,dy,dz,...
                      Hz_stagger_x,Hz_stagger_y,Hz_stagger_z,...
                      0,0,1,...
                      a_x,a_y,a_z,b_x,b_y,b_z);
    B1 = zeros(size(B3));
    B2 = zeros(size(B3));
    B = [B1,B2,B3];
end

function axb = crossProduct(a,b)
    c1 =  a(:,2).*b(:,3)-a(:,3).*b(:,2);
    c2 = -a(:,1).*b(:,3)+a(:,3).*b(:,1);
    c3 =  a(:,1).*b(:,2)-a(:,2).*b(:,1);
    axb = [c1,c2,c3];
end

function aTb = dotProduct(a,b)
    aTb = a(:,1).*b(:,1)+a(:,2).*b(:,2)+a(:,3).*b(:,3);
end

function x_new = pushParticle(x_old, v, dt)
    x_new = x_old + dt.*v;
end

function v_new = updateVelocity(x,v,q,m,E,B,kappa,dt)
    v_new = updateVelocityBoris(x,v,q,m,E,B,kappa,dt);
%     v_new = vayMethod(v,q,m,E,B,kappa,dt);
end

function u_next = updateVelocityBoris(x,u,q,m,E,B,kappa,dt)

    gamma_prev = computeGammaU(u,kappa);

    B_mag = sqrt(dotProduct(B,B));

    B_mag(B_mag == 0) = 1;

    B_hat = B ./ B_mag;

    u_minus = u + (q.*dt)./(2*m).*E;
    
    gamma_minus = computeGammaU(u_minus,kappa);
%     gamma_ave = (gamma_prev + gamma_next)/2;

    t = tan((q.*dt)./(2.*m.*gamma_minus).*B_mag).*B_hat;

    t_star = (q.*dt)./(2.*m.*gamma_minus).*B;

    t_mag = dotProduct(t,t);

    s = 2./(1 + t_mag).*t;
    
    u_prime = u_minus + crossProduct(u_minus,t);
    u_plus = u_minus + crossProduct(u_prime, s);
    u_next = u_plus + (q.*dt)./(2*m).*E;
% 
%     gam_minus = computeGammaU(u_minus,kappa);
%     gam_plus = computeGammaU(u_plus,kappa);
%     gam_next = computeGammaU(u_next,kappa);

%     assert(all(gam_minus - gam_plus < 1e-12));
%     assert(all(gam_plus - gam_next < 1e-12));

end

function u_next = vayMethod(u_prev,q,m,E,B,kappa,dt)
    
    gamma_prev = computeGammaU(u_prev,kappa);
    
    u_cross_B = crossProduct(u_prev,B);
    
    u_half = u_prev + (q.*dt)./(2.*m).*(E + u_cross_B./gamma_prev);
    
    u_prime = u_half + (q.*dt)./(2.*m).*E;
    
    tau = (q*dt)./(2.*m).*B;
    
    u_star = dotProduct(u_prime, tau./kappa);
    gamma_prime = computeGammaU(u_prime,kappa);
    sigma = gamma_prime.^2 - tau.^2;
    
    gamma_next = sqrt((sigma + sqrt(sigma.^2 + 4*(tau.^2 + u_star.^2)))./2);
    
    t = tau./gamma_next;
    
    t_mag2 = dotProduct(t,t);
    
    % Vay has 1/(t+t^2), where bold t is a vector. This isn't bolded though, so potential for error
    s = 1./(1 + t_mag2);
    
    u_next = s.*(u_prime + dotProduct(u_prime, t).*t + crossProduct(u_prime, t));

end

function gamma = computeGammaU(u,kappa)
    u_mag2 = dotProduct(u,u);
    gamma = sqrt(1+u_mag2./kappa^2);
end

function gamma = computeGammaV(v,kappa)
    v_mag2 = dotProduct(v,v);
    gamma = sqrt(1./(1-v_mag2./kappa^2));
end