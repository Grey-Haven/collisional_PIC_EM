clear;
q_ele = 1;
m_ele = 1;
tag_ele = '-';
p = Particle([0,0,0], [.1,0,0], tag_ele, q_ele, m_ele);

dt = .01;

a_x = -1;
a_y = -1;
a_z = -1;

b_x = 1;
b_y = 1;
b_z = 1;

zone_left = a_x;
zone_rite = b_x;
zone_bot = a_y;
zone_top = b_y;

L_x = b_x-a_x;
L_y = b_y-a_y;
L_z = b_z-a_z;

N = 64;

dx = (b_x-a_x)/(N);
dy = (b_y-a_y)/(N);
dz = dx;

x = a_x:dx:b_x;
y = a_y:dy:b_y;
z = a_z:dz:b_z;

Jx_stagger_x = dx/2;
Jx_stagger_y = 0;
Jx_stagger_z = 0;

Jy_stagger_x = 0;
Jy_stagger_y = dy/2;
Jy_stagger_z = 0;

Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;
Hz_stagger_z = 0;

Ex = zeros(N,N-1,N);
Ey = 0*ones(N-1,N,N);
Ez = 0*ones(N,N,N);

Jx = zeros(size(Ex));
Jy = zeros(size(Ey));

Hx = zeros(N,N,N);
Hy = zeros(N,N,N);
Hz = 1*ones(N,N,N);

particles = injectParticle(p);

for i = 1:1000
    boris_push;
%     [Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,0,L_y/2, ...
%                                    0,0,0,0, ...
%                                    0,0,0,0, ...
%                                    x,y,dx,dy,weightMap);
%     update_em_Hz_mode_dirichlet;
    loc = particles(1,1:3);
    scatter3(loc(1),loc(2),loc(3));
%     axis([-1 1 -1 1 -1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");
    hold on;
    drawnow;
end
hold off;
function particles = injectParticle(theParticle)

    particle_locations(1,:) = theParticle.location;
    particle_velocities(1,:) = theParticle.velocity;
    particle_tags(1) = double(theParticle.tag);
    particle_charges(1) = theParticle.charge;
    particle_masses(1) = theParticle.mass;
    particle_weights(1) = 1;

    particles = [particle_locations,particle_velocities,particle_tags,particle_masses,particle_weights,particle_charges];
end