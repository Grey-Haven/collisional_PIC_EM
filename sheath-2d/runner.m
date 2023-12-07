clear;
close all;
addpath(genpath([fileparts(pwd)]));
addpath(genpath([fileparts(pwd), '/common']));
addpath(genpath([fileparts(pwd), '/basic_boris']));
% profile -timestamp on

% switch over to Hz, Ex, Ey

% Nondimensional grid
grid_setup;

% Physical constants
M_electron = 9.109e-31; % [kg]
Q_electron = 1.602e-19; % [C] (intentionally positive)
c = 2.99792458e8;
mu_0 = 1.25663706e-6;
eps_0 = 8.854187817e-12;
k_B = 1.380649e-23; % Boltzmann constant

% Nondimensionalized units
M = M_electron; % [kg]
Q = Q_electron; % [C]


q_elec = -Q_electron/Q;
q_ion = Q_electron/Q;
m_elec = M_electron/M;
m_ion = 1836 * m_elec;

eleTag = '-';
ionTag = '+';
 
n_bar = 10^13; % In m^-3
% n_bar = 2.5*10^12;
T_e = 10000; % In Kelvin [K]

w_p = sqrt((n_bar*Q^2)/(M*eps_0)); % angular frequency
lambda_D = sqrt((eps_0*k_B*T_e)/(n_bar*Q^2)); % debye length

% L = lambda_D; % In meters [m]
n_13 = 10^13; % Originally the n_bar, lowered n_bar to increase sheath width, need to keep domain height/width the same
L = sqrt((eps_0*k_B*T_e)/(n_13*Q^2));


% Change domain to .05, trickle through everything else, can save cell count

% % More nondimensionalized units
T = 1/w_p;  % In seconds/radians [s/r]
V = L/T; % In [m/s] (thermal velocity lam_D*w_p)
kappa = c/V;

% nondimensionalization parameters for involutions
sig_1 = (M*eps_0)/(Q^2*T^2*n_bar);
sig_2 =  mu_0*Q^2*L^2*n_bar/M;

hydrogen_collisional_cross_section = .36e-9 / L; % in meters

dt = 1/(kappa*sqrt(1/dx^2+1/dy^2));

% N_particles = Nx*Ny*50; % want ~50 particles per cell
N_particles = 204800;
% N_particles = 25600;
% particle_count_scaling = 2;
% N_particles = particle_count_scaling*2.5e5;
w_original = ((b_x-a_x)*(b_y-a_y)*1)/N_particles;

T_final = 60;
% N_steps = 100001;
N_steps = ceil(T_final/dt);

Jx = zeros(size(Ex));
Jy = zeros(size(Ey));

Jx_stagger_x = dx/2;
Jx_stagger_y = 0;
Jy_stagger_x = 0;
Jy_stagger_y = dy/2;

Jx_padding_x = 1;
Jx_padding_y = 0;
Jy_padding_x = 0;
Jy_padding_y = 1;

Hz_stagger_x = dx/2;
Hz_stagger_y = dy/2;

Ex_x = a_x+dx/2:dx:b_x-dx/2;
Ex_y = y;
Ey_x = x;
Ey_y = a_y+dy/2:dy:b_y-dy/2;
Hz_x = a_x+dx/2:dx:b_x-dx/2;
Hz_y = a_y+dy/2:dy:b_y-dy/2;
    
%++++++++++++++++
% HALF STEP
%++++++++++++++++

[CHx, CHy] = computeCurlH_xy_no_boundaries(Hz, dx, dy);

% =====================================
% Update H fields from Curl of E fields
% =====================================
Ex = Ex + -dt/2*kappa^2*(CHx - sig_2*Jx);
Ey = Ey + -dt/2*kappa^2*(CHy - sig_2*Jy);
%+++++++++++++++++++++
%END HALF STEP SECTION
%+++++++++++++++++++++

w_elec = w_original;
w_ion = w_original;

electrons = inject_electrons(N_particles, eleTag, q_elec, m_elec, w_elec, ...
                             a_x, b_x, a_y, b_y);
ions = inject_ions(N_particles, ionTag, q_ion, m_ion, w_ion, ...
                             a_x, b_x, a_y, b_y);

col_type = 'collisionless';
col_title = 'Collisionless';

run = now;
projectRoot = "C:\Users\StWhite\Documents\stuff\MSU\advisor_work\repos\collisional_PIC_EM\sheath-2d";
vidPath = projectRoot + "\results\" + col_type + "\";
disp(vidPath);
vidName = "sheath_nondimen_debye" + Nx + "x" + Ny + ".mp4";
% csvPath = projectRoot + "\results\scatter_current_path\csv_files\" + run + "\";
vidObj = VideoWriter(vidPath + vidName, 'MPEG-4');
open(vidObj);

particles = [electrons; ions];
rho = scatter_rho(particles(:,1:3),particles(:,10),particles(:,9),x,y);
rho(1,:) = 0;
rho(end,:) = 0;
rho(:,1) = 0;
rho(:,end) = 0;

% A_inv = construct_poisson_solver(size(rho,1)-2);
% A = construct_poisson_matrix(size(rho,1)-2);

% phi_normal = -poisson_solver(rho(2:end-1,2:end-1)/sig_1,dx,dy,A_inv);
phi_normal = -poisson_solver(rho(2:end-1,2:end-1)/sig_1,dx,dy);

N = size(rho,1)-2;

assert(dx == dy);

rho_vec = reshape(rho(2:end-1,2:end-1),N^2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we need to use Gauss-Seidel
% Right now using sparse matrices
% GS for if we get real expensive
% -------------------------------- %

% GS_TOL = (dx*dy)^2;
% phi_GS = gauss_seidel(dx*dy*rho_vec/sig_1,GS_TOL);
% phi_GS_vec = reshape(phi_GS,N^2,1); % For testing
% Aphi_GS_vec = A*phi_GS_vec;
% assert(norm(Aphi_GS_vec + dx*dy*rho_vec/sig_1) < GS_TOL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% For testing purposes (if the grid isn't too big)
% Aphi = -rho/sigma, where A is the discrete Laplacian operator
% phi_normal_vec = reshape(phi_normal,N^2,1);
% Aphi_vec = A*phi_normal_vec;
% assert(norm(Aphi_vec + dx*dy*rho_vec/sig_1) < 1e-10);

phi = zeros(size(phi_normal)+2);
phi(2:end-1,2:end-1) = phi_normal;

phi_x = diff(phi,1,1)/dx;
phi_y = diff(phi,1,2)/dy;

Ex(:,:) = phi_x(:,2:end-1);
Ey(:,:) = phi_y(2:end-1,:);

record_stride = 100;

% N_steps = 500;
heat_hist = zeros(ceil(N_steps/record_stride),1);
num_eles_hist = zeros(ceil(N_steps/record_stride),1);
s_hist = zeros(ceil(N_steps/record_stride),1);

for s=0:N_steps-1
    t = s*dt;

    prev_particle_locs = particles(:,1:3);
    % Boris push
    electrons = particles(particles(:,7) == eleTag,:);
    ions = particles(particles(:,7) == ionTag,:);
    % collide_particles;
    boris_push;
    % collide_takizuka;
    particles = [ions;electrons];
    % Remove particles outside boundaries
    [particles,prev_particle_locs] = remove_particles(particles, prev_particle_locs, a_x, b_x, a_y, b_y);
    % Compute current
    [Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y, ...
                                   Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y, ...
                                   Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y, ...
                                   x,y,dx,dy,dt);

    % Yee(t)
    update_em_Hz_mode_dirichlet;

    if (mod(s,record_stride) == 0)
        visualize;
    end
end

ts = dt*(0:record_stride:N_steps-1);
heat_hist = heat_hist * (M_electron*V^2/k_B);

figure;
plot(ts,heat_hist);
xlabel("t");
ylabel("temp in kelvin");
title("Temperature " + col_title + ", " + Nx + "x" + Ny);
saveas(gcf,vidPath + "temperature_" + Nx + "x" + Ny + ".jpg");

figure;
plot(ts,num_eles_hist);
xlabel("t");
ylabel("n_e");
title("Number of Electrons " + col_title + ", " + Nx + "x" + Ny);
saveas(gcf,vidPath + "electron_count_" + Nx + "x" + Ny + ".jpg");

close(vidObj);
writematrix([ts',heat_hist],vidPath + 'heat_hist_' + Nx + "x" + Ny + '.csv');
writematrix([ts',num_eles_hist],vidPath + 'ele_hist_' + Nx + "x" + Ny + '.csv');


function writeCSVfiles(Ex,Ey,...
                       Jx,Jy,Bz,...
                       Ex_x,Ex_y,Ey_x,Ey_y,Bz_x,Bz_y,s,csv_path)
%     Ex_write = zeros(size(Ex)+1);
%     Ey_write = zeros(size(Ey)+1);
%     Bz_write = zeros(size(Bz)+1);
%     Jx_write = zeros(size(Jx)+1);
%     Jy_write = zeros(size(Jy)+1);
    
%     Ex_write(1,2:end) = Ex_x;
%     Ex_write(2:end,1) = Ex_y;
%     Ey_write(1,2:end) = Ey_x;
%     Ey_write(2:end,1) = Ey_y;
%     Bz_write(1,2:end) = Bz_x;
%     Bz_write(2:end,1) = Bz_y;
%     Jx_write(1,2:end) = Ex_x;
%     Jx_write(2:end,1) = Ex_y;
%     Jy_write(1,2:end) = Ey_x;
%     Jy_write(2:end,1) = Ey_y;
    
%     Ex_write(2:end,2:end) = Ex;
%     Ey_write(2:end,2:end) = Ey;
%     Bz_write(2:end,2:end) = Bz;
%     Jx_write(2:end,2:end) = Jx;
%     Jy_write(2:end,2:end) = Jy;
    
    Ex_path = csv_path + "Ex\";
    Ey_path = csv_path + "Ey\";
    Bz_path = csv_path + "Hz\";
    Jx_path = csv_path + "Jx\";
    Jy_path = csv_path + "Jy\";
    
    Ex_title = Ex_path + "Ex_" + num2str(s) + ".csv";
    Ey_title = Ey_path + "Ey_" + num2str(s) + ".csv";
    Bz_title = Bz_path + "Bz_" + num2str(s) + ".csv";
    Jx_title = Jx_path + "Jx_" + num2str(s) + ".csv";
    Jy_title = Jy_path + "Jy_" + num2str(s) + ".csv";    
    
    if ~isfolder(csv_path)
       mkdir(csv_path)
       mkdir(Ex_path);
       mkdir(Ey_path);
       mkdir(Bz_path);
       mkdir(Jx_path);
       mkdir(Jy_path);
    end

    Ex_table = matrix_to_table(Ex_x,Ex_y,Ex,'Ex');
    Ey_table = matrix_to_table(Ey_x,Ey_y,Ey,'Ey');
    Jx_table = matrix_to_table(Ex_x,Ex_y,Jx,'Jx');
    Jy_table = matrix_to_table(Ey_x,Ey_y,Jy,'Jy');
    Bz_table = matrix_to_table(Bz_x,Bz_y,Bz,'Bz');
    
    writetable(Ex_table,Ex_title);
    writetable(Ey_table,Ey_title);
    writetable(Bz_table,Bz_title);
    writetable(Jx_table,Jx_title);
    writetable(Jy_table,Jy_title);
    
end

function theTable = matrix_to_table(x,y,A,col3)
    T = zeros(size(A,1)*size(A,2),3);
    for i = 1:length(x)
        for j = 1:length(y)
            idx = (i-1)*size(A,1) + j;
            T(idx,1) = x(i);
            T(idx,2) = y(j);
            T(idx,3) = A(j,i);
        end
    end
    theTable = array2table(T,'VariableNames',{'x','y',col3});
end