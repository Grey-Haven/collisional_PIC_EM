clear;
M = 9.10938e-31;
Q = 1.602e-19;

sp_wt = 1;
mp_q = 1;

L_x = 4*pi;
L_y = 4*pi;

Nx = 32;
Ny = 32;

dx = L_x/Nx;
dy = L_y/Ny;

x = -L_x/2:dx:L_x/2-dx;
y = -L_y/2:dy:L_y/2-dx;

eleTag = '-';
w_elec = 1;

q_elec = -Q/Q;
m_elec = M/M;

weightMap = containers.Map({eleTag},{w_elec});
particles = Particle.empty(1,0);
theParticle = Particle([0,0,0],[.1,.1,0],eleTag,q_elec,m_elec);

Jx = zeros(Nx,Ny);
Jy = zeros(Nx,Ny);

Jx_stagger_x = 0;
Jx_stagger_y = 0;
Jy_stagger_x = 0;
Jy_stagger_y = 0;
% 
% Jx_padding_x = 0;
% Jx_padding_y = 0;
% Jy_padding_x = 0;
% Jy_padding_y = 0;
% 
% x0 = (x(15)+x(16))/2;
% y0 = (y(15)+y(16))/2;
% 
% particles{1} = theParticle;
% [Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,L_x/2,L_y/2,...
%                                Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
%                                Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
%                                x,y,dx,dy,weightMap);
% % surf(x,y,Jx);
% assert(size(Jx,1) == Nx && size(Jx,2) == Ny, "Jx size is incorrect");
% assert(size(Jy,1) == Nx && size(Jy,2) == Ny, "Jy size is incorrect");
% 
% assert(Jx(15,15) == Jx(15,16), "charge is not evenly distributed");
% assert(Jx(15,15) == Jx(16,15), "charge is not evenly distributed");
% assert(Jx(15,15) == Jx(16,16), "charge is not evenly distributed");
% assert(Jx(15,14) == 0, "should be zero elsewhere");
% 
% assert(Jy(15,15) == Jy(15,16), "charge is not evenly distributed");
% assert(Jy(15,15) == Jy(16,15), "charge is not evenly distributed");
% assert(Jy(15,15) == Jy(16,16), "charge is not evenly distributed");
% assert(Jy(15,14) == 0, "should be zero elsewhere");
% 
% 
% Jx = zeros(Nx,Ny);
% Jy = zeros(Nx,Ny);
% 
% x0 = x(15)*.1 + x(16)*.9;
% y0 = y(15)*.1 + y(16)*.9;
% theParticle.location = [x0, y0, 0];
% particles{1} = theParticle;
% [Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,L_x/2,L_y/2,...
%                                Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
%                                Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
%                                x,y,dx,dy,weightMap);
% assert(Jx(15,15) < Jx(15,16), "current is not unevenly distributed");
% assert(Jx(15,15) < Jx(16,15), "current is not unevenly distributed");
% assert(Jx(15,16) == Jx(16,15), "current is not unevenly distributed");
% assert(Jx(15,16) < Jx(16,16), "current is not unevenly distributed");
% 
% surf(x,y,Jx);
% 
% Jx = zeros(Nx,Ny);
% Jy = zeros(Nx,Ny);
% 
% x0 = x(2);
% y0 = y(2);
% theParticle.location = [x0, y0, 0];
% particles{1} = theParticle;
% 
% [Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,L_x/2,L_y/2,...
%                                Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
%                                Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
%                                x,y,dx,dy,weightMap);
% assert(Jx(2,2) == theParticle.charge*theParticle.velocity(1)*weightMap(eleTag)/(dx*dy), "current not correctly mapped");
% surf(x,y,Jx);
% 
% Jx = zeros(2,3);
% Jy = zeros(3,2);
% 
dx = .5;
dy = .5;

Jx_stagger_x = dx/2;
Jx_stagger_y = 0;
Jy_stagger_x = 0;
Jy_stagger_y = dy/2;

x = 0:dx:1;
y = -.5:dy:.5;
% 
% x0 = .51;
% y0 = .01;
% theParticle.location = [x0, y0, 0];
% particles{1} = theParticle;
% 
% [Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,0,.5,...
%                                Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
%                                Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
%                                x,y,dx,dy,weightMap);


Jx_padding_x = 1;
Jx_padding_y = 0;
Jy_padding_x = 0;
Jy_padding_y = 1;

Jx = zeros(2+Jx_padding_x*2,3+Jx_padding_y*2);
Jy = zeros(3+Jy_padding_x*2,2+Jy_padding_y*2);

x0 = .6;
y0 = -.1;
theParticle.location = [x0, y0, 0];
particles{1} = theParticle;

particles = injectParticle(theParticle);

[Jx_0,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,0,.5,...
                                 Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                                 Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                                 x,y,dx,dy);
x0 = .3;
y0 = -.1;
theParticle.location = [x0, y0, 0];
particles = injectParticle(theParticle);

[Jx_1,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,0,.5,...
                                 Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                                 Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                                 x,y,dx,dy);
assert(Jx_0(2,2) == Jx_1(2,2), "particles the same distance from a node should scatter the same value"); 

L_x = 1;
L_y = 1;

a_x = 0;
b_x = L_x;

a_y = -L_y/2;
b_y =  L_y/2;

a_z = 0;
b_z = 0;

Nx = 2;
Ny = 2;

dx = (b_x-a_x)/Nx;
dy = (b_y-a_y)/Ny;
dz = dx;

x = a_x:dx:b_x;
y = a_y:dy:b_y;
z = 0;

Hz_x_padding = 0;
Hz_y_padding = 0;

Ex_x_padding = 0;
Ex_y_padding = 0;

Ey_x_padding = 0;
Ey_y_padding = 0;

Ex_mag = .1;
Ey_mag = .2;
Hz_mag = 1;

% Hz mode
Hz = Hz_mag*ones(Nx+Hz_x_padding,Ny+Hz_y_padding);
Ex = Ex_mag*ones(Nx+Ex_x_padding,Ny-1+Ex_y_padding);
Ey = Ey_mag*ones(Nx-1+Ey_x_padding,Ny+Ey_y_padding);

% Ez mode
Hx = zeros(size(Ex));
Hy = zeros(size(Ey));
Ez = zeros(size(Hz));

Ex_padded = zeros(size(Ex,1)+2,size(Ex,2)+2,size(Ex,3));
Ex_padded(2:end-1,2:end-1,:) = Ex;

Ey_padded = zeros(size(Ey,1)+2,size(Ey,2)+2,size(Ey,3));
Ey_padded(2:end-1,2:end-1,:) = Ey;

Hz_padded = zeros(size(Hz,1)+2,size(Hz,2)+2,size(Hz,3));
Hz_padded(2:end-1,2:end-1,:) = Hz; % working for 2d, need to revisit for 3d probably

Ex_x = x;
Ex_y = y;
Ey_x = x;
Ey_y = y;

Ex_x_padded = Ex_x(1)-dx:dx:Ex_x(end)+dx;
Ey_y_padded = Ey_y(1)-dy:dy:Ey_y(end)+dy;
Hz_x_padded = x(1)-dx:dx:x(end)+dx;
Hz_y_padded = y(1)-dy:dy:y(end)+dy;

Ex_stagger_x = dx/2;
Ex_stagger_y = 0;
Ex_stagger_z = 0;
Ey_stagger_x = 0;
Ey_stagger_y = dy/2;
Ey_stagger_z = 0;
Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;
Hz_stagger_z = 0;

x0 = .05;
y0 = .1;
theParticle.location = [x0, y0, 0];
particles = injectParticle(theParticle);

[Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,0,.5,...
                               Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                               Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                               x,y,dx,dy);

disp(Jy);

Nx = 3;
Ny = 3;

dx = (b_x-a_x)/Nx;
dy = (b_y-a_y)/Ny;
dz = dx;

x = a_x:dx:b_x;
y = a_y:dy:b_y;
z = 0;

Jx = zeros(Nx,Ny-1);
Jy = zeros(Nx-1,Ny);

Jx_stagger_x = dx/2;
Jx_stagger_y = 0;
Jy_stagger_x = 0;
Jy_stagger_y = dy/2;

Hz_x_padding = 0;
Hz_y_padding = 0;

Ex_x_padding = 0;
Ex_y_padding = 0;

Ey_x_padding = 0;
Ey_y_padding = 0;

Ex_mag = .1;
Ey_mag = .2;
Hz_mag = 1;

% Hz mode
Hz = Hz_mag*ones(Nx+Hz_x_padding,Ny+Hz_y_padding);
Ex = Ex_mag*ones(Nx+Ex_x_padding,Ny-1+Ex_y_padding);
Ey = Ey_mag*ones(Nx-1+Ey_x_padding,Ny+Ey_y_padding);

% Ez mode
Hx = zeros(size(Ex));
Hy = zeros(size(Ey));
Ez = zeros(size(Hz));

Ex_padded = zeros(size(Ex,1)+2,size(Ex,2)+2,size(Ex,3));
Ex_padded(2:end-1,2:end-1,:) = Ex;

Ey_padded = zeros(size(Ey,1)+2,size(Ey,2)+2,size(Ey,3));
Ey_padded(2:end-1,2:end-1,:) = Ey;

Hz_padded = zeros(size(Hz,1)+2,size(Hz,2)+2,size(Hz,3));
Hz_padded(2:end-1,2:end-1,:) = Hz; % working for 2d, need to revisit for 3d probably

Ex_x = x(1:end-1) + dx/2;
Ex_y = y;
Ey_x = x;
Ey_y = y(1:end-1) + dy/2;

Ex_x_padded = Ex_x(1)-dx:dx:Ex_x(end)+dx;
Ey_y_padded = Ey_y(1)-dy:dy:Ey_y(end)+dy;
Hz_x_padded = x(1)-dx:dx:x(end)+dx;
Hz_y_padded = y(1)-dy:dy:y(end)+dy;

Ex_stagger_x = dx/2;
Ex_stagger_y = 0;
Ex_stagger_z = 0;
Ey_stagger_x = 0;
Ey_stagger_y = dy/2;
Ey_stagger_z = 0;
Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;
Hz_stagger_z = 0;

x0 = .1;
y0 = -.2;
theParticle.location = [x0, y0, 0];
particles = injectParticle(theParticle);

[Jx,Jy] = scatter_particles_xy(particles,Jx,Jy,L_x,L_y,0,.5,...
                               Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                               Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                               x,y,dx,dy);

Jy_with_boundaries = zeros(size(Jy,1)+2, size(Jy,2));
Jy_with_boundaries(2:end-1,:) = Jy;
Ey_y_with_boundaries = Ey_y(1)-dy:dy:Ey_y(end)+dy;

% charge*weight*velocity scattered to four adjacent nodes

d1 = (.1-0)/dx;
d2 = (-.2--(1/3))/dy;

q = theParticle.charge;
w = theParticle.mass;
vy = theParticle.velocity(2);

N1 = (1-d1)*(1-d2)*q*w*vy;
N2 = (d1)*(1-d2)*q*w*vy;
N3 = (1-d1)*(d2)*q*w*vy;
N4 = (d1)*(d2)*q*w*vy;

eps = 1e-16;

assert(Jy_with_boundaries(1,1) == 0); % PEC Boundaries
assert(abs(Jy_with_boundaries(2,1) - N2/(dx*dy)) < eps);
assert(Jy_with_boundaries(1,2) == 0); % PEC Boundaries
assert(abs(Jy_with_boundaries(2,2) - N4/(dx*dy)) < eps);

% surf(x(2:end-1),Ey_y_with_boundaries,Jy_with_boundaries');
% title("Jy");
% xlabel("x");
% ylabel("y");

function particles = injectParticle(theParticle)

    particle_locations(1,:) = theParticle.location;
    particle_velocities(1,:) = theParticle.velocity;
    particle_tags(1) = double(theParticle.tag);
    particle_charges(1) = theParticle.charge;
    particle_masses(1) = theParticle.mass;

    particles = [particle_locations,particle_velocities,particle_tags,particle_masses,particle_charges];
end