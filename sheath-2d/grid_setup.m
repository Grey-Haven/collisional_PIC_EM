L_x = 32;
L_y = 32;

offset_x = L_x/2;
offset_y = L_y/2;

a_x = -offset_x;
b_x = a_x + L_x;

a_y = -offset_y;
b_y = a_y + L_y;

a_z = 0;
b_z = 0;

Nx = 32;
Ny = 32;

dx = (b_x-a_x)/Nx;
dy = (b_y-a_y)/Ny;
dz = dx;

x = a_x:dx:b_x;
y = a_y:dy:b_y;
z = 0;

% Hz mode
Hz = zeros(Nx,Ny);
Ex = zeros(Nx,Ny-1);
Ey = zeros(Nx-1,Ny);

% Ez mode
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
Ez = zeros(Nx,Ny);

rho = zeros(Nx+1,Ny+1);
rng(2);