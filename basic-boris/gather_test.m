clear;
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

Ex_padded = zeros(size(Ex,1)+2,size(Ex,2)+2,size(Ex,3)+2);
Ex_padded(2:end-1,2:end-1,2:end-1) = Ex;

Ey_padded = zeros(size(Ey,1)+2,size(Ey,2)+2,size(Ey,3)+2);
Ey_padded(2:end-1,2:end-1,2:end-1) = Ey;

Hz_padded = zeros(size(Hz,1)+2,size(Hz,2)+2,size(Hz,3)+2);
Hz_padded(2:end-1,2:end-1,2:end-1) = Hz; % working for 2d, need to revisit for 3d probably

Ex_x = x;
Ex_y = y;
Ey_x = x;
Ey_y = y;

x_padded = x(1)-dx:dx:x(end)+dx;
y_padded = y(1)-dy:dy:y(end)+dy;
z_padded = z(1)-dz:dz:z(end)+dz;

Ex_stagger_x = dx/2;
Ex_stagger_y = 0;
Ex_stagger_z = 0;
Ey_stagger_x = 0;
Ey_stagger_y = dy/2;
Ey_stagger_z = 0;
Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;
Hz_stagger_z = 0;

Ex_x_padded = (Ex_x(1)-dx:dx:Ex_x(end)+dx)+Ex_stagger_x;
Ey_y_padded = (Ey_y(1)-dy:dy:Ey_y(end)+dy)+Ey_stagger_y;
Hz_x_padded = (x(1)-dx:dx:x(end)+dx)+Hz_stagger_x;
Hz_y_padded = (y(1)-dy:dy:y(end)+dy)+Hz_stagger_y;

loc = [.01, 0, 0];

E1 = gather_field(Ex_padded,loc,Ex_x_padded,y,z_padded,dx,dy,dz,...
                  Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);
% assert(E1 == Ex_mag)

loc = [x(end)-.01, 0, 0];
E1 = gather_field(Ex_padded,loc,Ex_x_padded,y,z_padded,dx,dy,dz,...
                  Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);

loc = [.5,0,0];

E2 = gather_field(Ey_padded,loc,x,Ey_y_padded,z_padded,dx,dy,dz,...
                  Ey_stagger_x,Ey_stagger_y,Ey_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(E2 == Ey_mag);

loc = [.01,.2,0];

E1 = gather_field(Ex_padded,loc,Ex_x_padded,y,z_padded,dx,dy,dz,...
                  Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(E1 == Ex_mag*.312);

loc = [.95,-.4,0];

E1 = gather_field(Ex_padded,loc,Ex_x_padded,y,z_padded,dx,dy,dz,...
                  Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(E1 == Ex_mag*.12);

loc = [.05,.1,0];
B3 = gather_field(Hz_padded,loc,Hz_x_padded,Hz_y_padded,z_padded,dx,dy,dz,...
                  Hz_stagger_x,Hz_stagger_y,Hz_stagger_z,...
                  0,0,1,...
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(B3 == Hz_mag*(.18+.42));

% TEST CASE FOR 32x32 (breaking in runner code)
Nx = 32;
Ny = 32;

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

Ex_padded = zeros(size(Ex,1)+2,size(Ex,2)+2,size(Ex,3)+2);
Ex_padded(2:end-1,2:end-1,2:end-1) = Ex;

Ey_padded = zeros(size(Ey,1)+2,size(Ey,2)+2,size(Ey,3)+2);
Ey_padded(2:end-1,2:end-1,2:end-1) = Ey;

Hz_padded = zeros(size(Hz,1)+2,size(Hz,2)+2,size(Hz,3)+2);
Hz_padded(2:end-1,2:end-1,2:end-1) = Hz;

Ex_stagger_x = dx/2;
Ex_stagger_y = 0;
Ex_stagger_z = 0;
Ey_stagger_x = 0;
Ey_stagger_y = dy/2;
Ey_stagger_z = 0;
Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;
Hz_stagger_z = 0;

Ex_x_padded = (Ex_x(1)-dx:dx:Ex_x(end)+dx)+Ex_stagger_x;
Ey_y_padded = (Ey_y(1)-dy:dy:Ey_y(end)+dy)+Ey_stagger_y;
Hz_x_padded = (x(1)-dx:dx:x(end)+dx)+Hz_stagger_x;
Hz_y_padded = (y(1)-dy:dy:y(end)+dy)+Hz_stagger_y;

Hz_x_padded = (Ex_x(1)-dx:dx:Ex_x(end)-dx)+dx/2;
Hz_y_padded = (Ey_y(1)-dy:dy:Ey_y(end)-dy)+dy/2;

Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;
Hz_stagger_z = 0;

loc = [0.9378, 0.0396, 0];
B3 = gather_field(Hz_padded,loc,Hz_x_padded,Hz_y_padded,z_padded,dx,dy,dz,...
                  Hz_stagger_x,Hz_stagger_y,Hz_stagger_z,...
                  0,0,1,...
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(B3 == Hz_mag);
E1 = gather_field(Ex_padded,loc,Ex_x_padded,y,z_padded,dx,dy,dz,...
                  Ex_stagger_x,Ex_stagger_y,Ex_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(E1 == Ex_mag);
E2 = gather_field(Ey_padded,loc,x,Ey_y_padded,z_padded,dx,dy,dz,...
                  Ey_stagger_x,Ey_stagger_y,Ey_stagger_z,...
                  0,0,1,... % padding offset on x,y,z, respectively
                  a_x,a_y,a_z,b_x,b_y,b_z);
assert(abs(E2 - Ey_mag) < eps);