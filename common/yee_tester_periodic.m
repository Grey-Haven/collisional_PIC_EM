clear;
addpath(genpath([fileparts(pwd)]));
addpath(genpath([fileparts(pwd), '/common']));
addpath(genpath([fileparts(pwd), '/basic_boris']));

c = 2.9979e8;
M_electron = 9.109e-31; % [kg]
Q_electron = 1.602e-19; % [C] (intentionally positive)

L = 1e-1; % In meters [m]
% V = 5e7;  % In meters per second [m/s]
V = c; % Getting an analytic solution requires kappa = 1 and alpha == beta
T = L/V;  % In seconds [s] (time for particle to cross the grid)

% Nondimensional grid
grid_setup;

v_injection = 1; % rate of injection
num_particle_crossings = 1000; % number of particle crossings

T_final = num_particle_crossings; % how long it takes for all the particles to cross

% subject to change, make sure dt <= CFL
N_steps = 1e5;
dt = T_final/N_steps;

mu_0 = 1.25663706e-6;
eps_0 = 8.8541878128e-12;

kappa = c/V; % Nondimensional wave speed

M = M_electron; % [kg]
Q = Q_electron; % [C]

q_elec = -Q_electron/Q;
m_elec = 1;

% number_density = 7.8025e14; % number density [m^3]
number_density = 1.5522581e+14; % smaller density [m^3]
N_particles_est = 6e4; % estimated particles in simulation

% nondimensionalization parameter
sig_1 = (M*eps_0)/(Q^2*T^2*number_density);
sig_2 =  mu_0*Q^2*L^2*number_density/M;

eleTag = '-';
dt = 1/(2*kappa*sqrt(1/dx^2+1/dy^2));

injection_setup;

steps_per_crossing =  L_x/(v_injection*dt);
N_part_max = injection_rate*steps_per_crossing; % rough max of particles at any given time
w_original = (L_x*zone_height*1)/N_part_max;

ramp_start_time = 0;
ramp_end_time = steps_per_crossing*dt;

iterations = 5;

final_times = zeros(1,iterations+1);
diffs = zeros(1,iterations+1);
NX = 8;
NY = 8;

for iterate=0:iterations
    
    d_mult = 2^(iterate);

    disp(strcat('d: ', num2str(d_mult)));

    Nx = d_mult*NX;
    Ny = d_mult*NY;

    a_x = 0;
    b_x = L_x;

    a_y = 0;
    b_y =  L_y;

    a_z = 0;
    b_z = 0;

    dx = (b_x-a_x)/Nx;
    dy = (b_y-a_y)/Ny;
    dz = dx;

    x = a_x:dx:b_x - dx;
    y = a_y:dy:b_y - dy;
    z = 0;

    % Ez mode
    % Hx = zeros(Nx+1,Ny);
    % Hy = zeros(Nx,Ny+1);
    % Ez = zeros(Nx,Ny);

    % Hz mode
    Hz = zeros(Nx,Ny);
    Ex = zeros(Nx,Ny);
    Ey = zeros(Nx,Ny);

    Jx = zeros(size(Ex));
    Jy = zeros(size(Ey));

    Jx_stagger_x = dx/2;
    Jx_stagger_y = 0;
    Jy_stagger_x = 0;
    Jy_stagger_y = dy/2;

    Jx_padding_x = 0;
    Jx_padding_y = 0;
    Jy_padding_x = 0;
    Jy_padding_y = 0;

    Hz_stagger_x = dx/2;
    Hz_stagger_y = dy/2;

    Ex_x = a_x+dx/2:dx:b_x-dx/2;
    Ex_y = y;
    Ey_x = x;
    Ey_y = a_y+dy/2:dy:b_y-dy/2;
    Hz_x = Ex_x;
    Hz_y = Ey_y;
    
    STEPS = 250*d_mult;
    dt = 1/(2*kappa*sqrt(1/dx^2+1/dy^2));
    tau = dt*10;
    t_0 = 6*tau;
    x_src_idx = 10*d_mult;
    y_src_idx = 10*d_mult;
    
    n = 1;
    m = 1;
    E_0 = 1;
    H_0 = 1;
    alpha = 2*n*pi/b_x;
%     beta =  2*m*pi/b_y;
    beta = alpha;
    omega = sqrt(alpha^2+beta^2);
    [Ex,Ey,Hz] = initial_conditions(E_0,alpha,beta,omega,kappa,...
                                    Ex_x,Ex_y,...
                                    Ey_x,Ey_y,...
                                    Hz_x,Hz_y,dt);

%     Ex = zeros(size(Ex));
%     Ey = zeros(size(Ey));
%     Hz = zeros(size(Hz));
                                
    for s = 1:STEPS

        t = s*dt;
        update_em_Hz_mode_periodic;
%         A = analytic_solution_Ex(H_0,alpha,beta,omega,kappa,Hz,dx,dy,t);
        if (mod(s,10) == 0)
            Ex_A = analytic_solution_Ex(H_0,alpha,beta,omega,kappa,Ex_x,Ex_y,t);
            Ey_A = analytic_solution_Ey(H_0,alpha,beta,omega,kappa,Ey_x,Ey_y,t);
            Hz_A = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz_x,Hz_y,t);

            Ex_with_boundaries = zeros(size(Ex,1),size(Ex,2)+2);
            Ex_with_boundaries(:,2:end-1) = Ex(:,:);
            Ex_with_boundaries(:,1) = analytic_solution_Ex(H_0,alpha,beta,omega,kappa,Ex_x,y(1),t);
            Ex_with_boundaries(:,end) = analytic_solution_Ex(H_0,alpha,beta,omega,kappa,Ex_x,y(end),t);
            
            Ey_with_boundaries = zeros(size(Ey,1)+2,size(Ey,2));
            Ey_with_boundaries(2:end-1,:) = Ey(:,:);
            Ey_with_boundaries(1,:) = analytic_solution_Ey(H_0,alpha,beta,omega,kappa,x(1),Ey_y,t);
            Ey_with_boundaries(end,:) = analytic_solution_Ey(H_0,alpha,beta,omega,kappa,x(end),Ey_y,t);

%             subplot(2,1,1);
%             surf(Ex_x,Ex_y,Ex_with_boundaries');
%             title("Ex Numeric");
%             zlim([-1,1]);
%             xlabel("x");
%             ylabel("y");
%             caxis([-1,1]);
%             colorbar;

%             subplot(2,1,1);
%             surf(Ey_x,Ey_y,Ey');
%             title("Ey Numeric");
%             zlim([-1,1]);
%             xlabel("x");
%             ylabel("y");
%             caxis([-1,1]);
%             colorbar;

            subplot(2,1,1);
            surf(Hz_x,Hz_y,Hz');
            title("Hz Numeric");
            zlim([-1,1]);
            xlabel("x");
            ylabel("y");
            caxis([-1,1]);
            colorbar;
%             
%             subplot(2,1,2);
%             surf(Ex_x,Ex_y,Ex_A');
%             title("Ex Analytic");
%             zlim([-1,1]);
%             xlabel("x");
%             ylabel("y");
%             caxis([-1,1]);
%             colorbar;
%             
%             subplot(2,1,2);
%             surf(Ey_x,Ey_y,Ey_A');
%             title("Ey Analytic");
%             zlim([-1,1]);
%             xlabel("x");
%             ylabel("y");
%             caxis([-1,1]);
%             colorbar;
% 
            subplot(2,1,2);
            surf(Hz_x,Hz_y,Hz_A');
            title("Hz Analytic");
            zlim([-1,1]);
            xlabel("x");
            ylabel("y");
            caxis([-1,1]);
            colorbar;
        sgtitle("STEP: " + num2str(s) + ", t = " + num2str(t));
            drawnow;
        end

    end
    Hz_A = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz_x,Hz_y,t);
    final_times(iterate+1) = t;
    dts(iterate+1) = dt;
    for i = 1:size(Hz_A,1)
        for j = 1:size(Hz_A,2)
            diffs(iterate+1) = diffs(iterate+1) + abs(Hz_A(i,j) - Hz(i,j))^2;
        end
    end
    diffs(iterate+1) = sqrt(dx*dy*diffs(iterate+1));
end
figure;
loglog(dts,diffs);
linspace = dts;
linear = linspace;
quadratic = linspace.^2;
cubic = linspace.^3;
hold on;
loglog(linspace,linear);
loglog(linspace,quadratic);
loglog(linspace,cubic);
legend("Errors", "linear", "quadratic", "cubic");
title("Difference between analytic and numeric solution"); 
disp("Final Times:");
disp(final_times);
hold off;
%close(vidObj);

% ========================
% End Program
% ========================
function g_E = energy_source(t, t_0, tau)
    g_E = exp(-((t-t_0)/tau)^2);
end


function [Ex,Ey,Hz] = initial_conditions(H_0,alpha,beta,omega,kappa,Ex_x,Ex_y,Ey_x,Ey_y,Hz_x,Hz_y,dt)

    % =======================================
    % Initialize Ex and Ey field at t = dt/2
    % Initialize Hz field at t = 0
    % =======================================
    Hz = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz_x,Hz_y,0);
    Ex = analytic_solution_Ex(H_0,alpha,beta,omega,kappa,Ex_x,Ex_y,dt/2);
    Ey = analytic_solution_Ey(H_0,alpha,beta,omega,kappa,Ey_x,Ey_y,dt/2);
end

function analytic = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,x,y,t)
    analytic = zeros(length(x),length(y));
    for i = 1:length(x)
        xi = x(i);
        for j = 1:length(y)
            yj = y(j);
            analytic(i,j) = -(2*H_0/(kappa*sqrt(2)))*sin(alpha*xi)*sin(beta*yj)*cos(omega*kappa*t);
        end
    end
end

function analytic = analytic_solution_Ex(E_0,alpha,beta,omega,kappa,x,y,t)
    analytic = zeros(length(x),length(y));
    for i = 1:length(x)
        for j = 1:length(y)
            xi = x(i);
            yj = y(j);
            analytic(i,j) = -E_0*sin(alpha*xi)*cos(beta*yj)*sin(kappa*omega*t);
        end
    end
end

function analytic = analytic_solution_Ey(E_0,alpha,beta,omega,kappa,x,y,t)
    analytic = zeros(length(x),length(y));
    for i = 1:length(x)
        for j = 1:length(y)
            xi = x(i);
            yj = y(j);
            analytic(i,j) = E_0*cos(alpha*xi)*sin(beta*yj)*sin(kappa*omega*t);
        end
    end
end