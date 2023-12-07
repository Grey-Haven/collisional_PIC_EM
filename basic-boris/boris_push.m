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

particle_locations = electrons(:,1:3);
particle_velocities = electrons(:,4:6);
particle_masses = electrons(:,8);
particle_charges = electrons(:,10);

particles_in_domain_idxs = find(particle_in_zone(electrons(:,:), a_x, b_x, a_y, b_y));
particles_out_domain_idxs = find(~particle_in_zone(electrons(:,:), a_x, b_x, a_y, b_y));

particles_in_domain = electrons(particles_in_domain_idxs,:);
particles_out_domain = electrons(particles_out_domain_idxs,:);

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
    newVelocity = updateVelocity(xe,ve,qe,me,E,B,dt);
    particles_in_domain(:,4:6) = newVelocity;
    newLocation = pushParticle(xe,newVelocity,dt);
    particles_in_domain(:,1:3) = newLocation;
end
if (size(particles_out_domain,1) ~= 0)

    xe = particles_out_domain(:,1:3);
    ve = particles_out_domain(:,4:6);
    me = particles_out_domain(:,8);
    qe = particles_out_domain(:,10);

    newLocation = pushParticle(xe,ve,dt);
    particles_out_domain(:,1:3) = newLocation;
end

electrons = [particles_in_domain; particles_out_domain];

% for i = 1:size(particle_locations,1)-injection_rate
%     xe = particle_locations(i,:);
%     ve = particle_velocities(i,:);
%     qe = particle_charges(i);
%     me = particle_masses(i);
%     E = evalE(Ex_padded,Ey_padded,Ez_padded,xe,...
%               x,y,z,x_padded,y_padded,z_padded,...
%               dx,dy,dz,...
%               Jx_stagger_x,Jx_stagger_y,0,...
%               Jy_stagger_x,Jy_stagger_y,0,...
%               0,0,0,... % no Ez field, no staggering
%               a_x,a_y,a_z,b_x,b_y,b_z);
%     B = evalB(Hx_padded,Hy_padded,Hz_padded,xe,...
%               x,y,z,x_padded,y_padded,z_padded,...
%               dx,dy,dz,...
%               0,0,0,... % no Hx field
%               0,0,0,... % no Hy field
%               Hz_stagger_x,Hz_stagger_y,0,...
%               a_x,a_y,a_z,b_x,b_y,b_z);
%     newVelocity = updateVelocity(xe,ve,qe,me,E,B,dt);
%     particle_velocities(i,:) = newVelocity(:);
%     newLocation = pushParticle(xe,ve,dt);
% %     if ~(abs(newLocation(2)-particle_locations(i,2)) < dy)
% %         disp('foo');
% %     end
% %     if ~(abs(newLocation(1)-particle_locations(i,1)) < dx)
% %         disp('foo');
% %     end
%     assert(abs(newLocation(1)-particle_locations(i,1)) < dx);
%     assert(abs(newLocation(2)-particle_locations(i,2)) < dy);
%     particle_locations(i,:) = newLocation(:);
% end
% for i = size(particle_locations,1)-injection_rate+1:size(particle_locations,1)
%     xe = particle_locations(i,:);
%     ve = particle_velocities(i,:);
%     qe = particle_charges(i);
%     me = particle_masses(i);
%     newLocation = pushParticle(xe,ve,dt);
%     assert(abs(newLocation(1)-particle_locations(i,1)) < dx);
%     assert(abs(newLocation(2)-particle_locations(i,2)) < dy);
%     particle_locations(i,:) = newLocation(:);
% end
% particles(:,1:3) = particle_locations;
% particles(:,4:6) = particle_velocities;


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

function x_new = pushParticle(x_old, v, dt)
    x_new = x_old + v*dt;
end

function v_new = updateVelocity(x,v,q,m,E,B,dt)
    v_new = updateVelocityBoris(x,v,q,m,E,B,dt);
end

function v_b = updateVelocityBoris(x,v,q,m,E,B,dt)
    v_b = zeros(size(v));
    v_minus = zeros(size(v));
    v_plus  = zeros(size(v));
    v_prime = zeros(size(v));
    
    t = zeros(size(v));
    s = zeros(size(v));
    
    for dim=1:3
        t(:,dim) = q./m.*B(:,dim)*.5*dt;
    end
    
    t_mag2 = t(:,1).*t(:,1)+t(:,2).*t(:,2)+t(:,3).*t(:,3);
    
    for dim=1:3
        s(:,dim) = 2*t(:,dim)./(1+t_mag2);
    end
    
    for dim=1:3
        v_minus(:,dim) = v(:,dim) + q./m.*E(:,dim)*.5*dt;
    end
    
    v_minus_cross_t = crossProduct(v_minus, t);
    
    for dim=1:3
        v_prime(:,dim) = v_minus(:,dim) + v_minus_cross_t(:,dim);
    end
    
    v_prime_cross_s = crossProduct(v_prime, s);
    
    for dim=1:3
        v_plus(:,dim) = v_minus(:,dim) + v_prime_cross_s(:,dim);
    end
    
    for dim=1:3
        v_b(:,dim) = v_plus(:,dim) + q./m.*E(:,dim)*.5*dt;
    end
    
end