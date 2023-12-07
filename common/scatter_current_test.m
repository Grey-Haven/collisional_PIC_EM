clear;
r = [0.05,.1,0];
v = [1,0,0];
eleTag = 45;
q_elec = -1;
m_elec = 1;
w_elec = 2.9e-5; % roughly what w_original is

L_x = 1;
L_y = 1;

Nx = 2;
Ny = 2;

dx = L_x/Nx;
dy = L_y/Ny;

x = 0:dx:L_x;
y = -L_y/2:dy:L_y/2;

dt = .0018;

Jx = zeros(Nx,Ny-1);
Jy = zeros(Nx-1,Ny);

Jx_stagger_x = dx/2;
Jx_stagger_y = 0;
Jy_stagger_x = 0;
Jy_stagger_y = dy/2;

Jx_padding_x = 1;
Jx_padding_y = 0;
Jy_padding_x = 0;
Jy_padding_y = 1;

particles = [r,v,eleTag,m_elec,w_elec,q_elec];
prev_particle_locs = [-0.05,.1,0];

offset_y = L_y/2;
offset_x = 0;

[Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
                               Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                               Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                               x,y,dx,dy,dt);

                           
r = [.45,.1,0];
% d = vt
d = r-prev_particle_locs;
v = d/dt;
particles = [r,v,eleTag,m_elec,w_elec,q_elec];
prev_particle_locs = [.4,.1,0];

J = q_elec*v;

[Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
                               Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                               Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                               x,y,dx,dy,dt);
         
disp(Jx);                  
disp(Jy);
disp(J);
                           
                           
                           
% r = [.6,.1,0];
% particles = [r,v,eleTag,m_elec,w_elec,q_elec];
% prev_particle_locs = [.4,-.1,0]; % Tougher edge case, fix later
% 
% [Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
%                                Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
%                                Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
%                                x,y,dx,dy,dt);
                           
                           
                           
r = [.6,.2,0];
particles = [r,v,eleTag,m_elec,w_elec,q_elec];
prev_particle_locs = [.4,-.1,0];

offset_y = L_y/2;
offset_x = 0;

[Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
                               Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                               Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                               x,y,dx,dy,dt);

Nx = 4;
Ny = 4;

dx = L_x/(Nx);
dy = L_y/(Ny);

x = 0:dx:L_x;
y = -L_y/2:dy:L_y/2;

dt = .0018;

Jx = zeros(Nx,Ny-1);
Jy = zeros(Nx-1,Ny);

% Staying in the same cell
r = [.4,-.1,0];
particles = [r,v,eleTag,m_elec,w_elec,q_elec];
prev_particle_locs = [.3,-.2,0];

% Jx = zeros(size(Jx,1)+Jx_padding_x*2, size(Jx,2)+Jx_padding_y*2);
% Jy = zeros(size(Jy,1)+Jy_padding_x*2, size(Jy,2)+Jy_padding_y*2);

[Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
                               Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                               Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                               x,y,dx,dy,dt);
