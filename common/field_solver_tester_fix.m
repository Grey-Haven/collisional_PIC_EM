clear;
addpath(genpath([fileparts(pwd)]));
addpath(genpath([fileparts(pwd), '/common']));
addpath(genpath([fileparts(pwd), '/basic_boris']));

c = 2.9979e8;
M_electron = 9.109e-31; % [kg]
Q_electron = 1.602e-19; % [C] (intentionally positive)

L = 1e-1; % In meters [m]
% V = 5e7;  % In meters per second [m/s]
V = c;
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
dt = 1/(kappa*sqrt(1/dx^2+1/dy^2));

injection_setup;

steps_per_crossing =  L_x/(v_injection*dt);
N_part_max = injection_rate*steps_per_crossing; % rough max of particles at any given time
w_original = (L_x*zone_height*1)/N_part_max;

ramp_start_time = 0;
ramp_end_time = steps_per_crossing*dt;

iterations = 5;

final_times = zeros(1,iterations+1);

NX = 64;
NY = 64;

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

    dx = (b_x-a_x)/(Nx);
    dy = (b_y-a_y)/(Ny);
    dz = dx;

    x = a_x:dx:b_x-dx;
    y = a_y:dy:b_y-dy;
    z = 0;

    % Ez mode
    % Hx = zeros(Nx+1,Ny);
    % Hy = zeros(Nx,Ny+1);
    % Ez = zeros(Nx,Ny);

    Hz_x_padding = 0;
    Hz_y_padding = 0;

    Ex_x_padding = 0;
    Ex_y_padding = 0;

    Ey_x_padding = 0;
    Ey_y_padding = 0;

    % Hz mode
    Hz = zeros(Nx+Hz_x_padding,Ny+Hz_y_padding);
    Ex = zeros(Nx+Ex_x_padding,Ny-1+Ex_y_padding);
    Ey = zeros(Nx-1+Ey_x_padding,Ny+Ey_y_padding);

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

    Ex_x = a_x-dx/2:dx:b_x+dx/2;
    Ey_y = a_y-dy/2:dy:b_y+dy/2;
    
    STEPS = 50*d_mult;
    dt = 1/(2*kappa*sqrt(1/dx^2+1/dy^2));
    tau = dt*10;
    t_0 = 6*tau;
    x_src_idx = 10*d_mult;
    y_src_idx = 10*d_mult;

    n = 1;
    m = 1;
    E_0 = 1;
    H_0 = 1;
    alpha = n*pi/b_x;
    beta =  m*pi/b_y;
    omega = sqrt(alpha^2+beta^2);

    [Ex,Ey,Hz] = initial_conditions(E_0,alpha,beta,omega,kappa,dx,dy,dt,Ex,Ey,Hz);
%     A = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz,dx,dy,0);
%     subplot(2,1,1);
%     surf(Ex_x,Ey_y,Hz');
%     subplot(2,1,2);
%     surf(Ex_x,Ey_y,A');
        
%     Ey(1,:) = 0;
%     Ey(end,:) = 0;
%     Ey(:,1) = 0;
%     Ey(:,end) = 0;
%     Ex(1,:) = 0;
%     Ex(end,:) = 0;
%     Ex(:,1) = 0;
%     Ex(:,end) = 0;
    
    for s = 1:STEPS

        
%         Hz(2,:) = 0;
%         Hz(end-1,:) = 0;
%         Hz(:,2) = 0;
%         Hz(:,end-1) = 0;
%         
%         Ey(1,:) = 0;
%         Ey(end,:) = 0;
%         Ey(:,1) = 0;
%         Ey(:,end) = 0;
%         Ex(1,:) = 0;
%         Ex(end,:) = 0;
%         Ex(:,1) = 0;
%         Ex(:,end) = 0;

        t = s*dt;
        update_em_Hz_mode_dirichlet_fix;
%         Hz(1,:) = 0;
%         Hz(Nx,:) = 0;
%         Hz(:,1) = 0;
%         Hz(:,Ny) = 0;
%         Ex(:,1) = 0;
%         Ex(:,end) = 0;
%         Ey(1,:) = 0;
%         Ey(end,:) = 0;
%         surf(x,y,Hz');
%         zlim([-1,1]);
%         drawnow;
%         A = analytic_solution_Ey(H_0,alpha,beta,omega,kappa,Hz,dx,dy,t);
        if (mod(s,10) == 0)
            A = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz,dx,dy,t);
            subplot(2,1,1);
            surf(x,y,Hz');
            title("Hz Numeric");
            zlim([-1,1]);
            xlabel("x");
            ylabel("y");
            caxis([-1,1]);
            colorbar;
            
            subplot(2,1,2);
            surf(x,y,A');
            title("Hz Analytic");
            zlim([-1,1]);
            xlabel("x");
            ylabel("y");
            caxis([-1,1]);
            colorbar;
            drawnow;
        end
%         Hz(x_src_idx,y_src_idx) = Hz(x_src_idx,y_src_idx) + energy_source(t, t_0, tau);
%         if (mod(s,500) == 0)
%             A = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz,dx,dy,t);
%             surf(Ex_x,Ey_y,A');
%             shading interp;
%             title("Hz");
%             xlim([Ex_x(1),Ex_x(end)]);
%             ylim([Ey_y(1),Ey_y(end)]);
%             colorbar;
%             caxis([-.01,.01]);
%             view(2);
%             drawnow;
%         end

    end
    
    if (d_mult ~= 1)
        curr_results = Hz(2:2^(iterate):end,2:2^(iterate):end);
        diff = curr_results-prev_results;
        sum_diff = 0;
        max_diff = 0;
        for ix = 1:size(curr_results,1)
            for jy = 1:size(curr_results,2)
                sum_diff = sum_diff + diff(ix,jy)^2;
                max_diff = max(max_diff,abs(diff(ix,jy)));
            end
        end
        %diffs(iterate) = sqrt(sum_diff);
        diffs(iterate) = max_diff;
        prev_results=curr_results;
        dts(iterate) = dt;
        final_times(iterate+1) = t;
    else
        prev_results = Hz(1:end,1:end);
        final_times(iterate+1) = t;
    end
    
end
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
title("Difference between dt and previous dt"); 
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


function [Ex,Ey,Hz] = initial_conditions(H_0,alpha,beta,omega,kappa,dx,dy,dt,Ex,Ey,Hz)

    % =======================================
    % Initialize Ex and Ey field at t = dt/2
    % Initialize Hz field at t = 0
    % =======================================
    Hz = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz,dx,dy,0);
    Ex = analytic_solution_Ex(H_0,alpha,beta,omega,kappa,Ex,dx,dy,dt/2);
    Ey = analytic_solution_Ey(H_0,alpha,beta,omega,kappa,Ey,dx,dy,dt/2);
end

function analytic = analytic_solution_Hz(H_0,alpha,beta,omega,kappa,Hz,dx,dy,t)
    analytic = zeros(size(Hz));
    for i = 1:size(Hz,1)
        xi_stagger = (i-1)*dx + dx/2;
        for j = 1:size(Hz,2)
            yj_stagger = (j-1)*dy + dy/2;
%             analytic(i,j) = H_0*sin(alpha*xi_stagger)*sin(beta*yj_stagger)*cos(kappa*omega*t);
            analytic(i,j) = H_0*cos(alpha*xi_stagger)*cos(beta*yj_stagger)*cos(omega*kappa*t);
        end
    end
end

function analytic = analytic_solution_Ex(E_0,alpha,beta,omega,kappa,Ex,dx,dy,t)
    analytic = zeros(size(Ex));
    for i = 1:size(Ex,1) % padding
        for j = 1:size(Ex,2)
            xi = (i-1)*dx;
            yj = (j)*dy;
            xi_stagger = xi + dx/2;
            yj_stagger = yj + dy/2;
%             analytic(i,j) = beta*kappa/omega*E_0*sin(alpha*xi_stagger)*cos(beta*yj)*sin(kappa*omega*t);
            analytic(i,j) = E_0*cos(alpha*xi_stagger)*sin(beta*yj)*sin(kappa*omega*t);
        end
    end
%     analytic(:,1) = 0;
%     analytic(:,end) = 0;
end

function analytic = analytic_solution_Ey(E_0,alpha,beta,omega,kappa,Ey,dx,dy,t)
    analytic = zeros(size(Ey));
    for i = 1:size(Ey,1)
        for j = 1:size(Ey,2) % padding
            xi = (i)*dx;
            yj = (j-1)*dy;
            xi_stagger = xi + dx/2;
            yj_stagger = yj + dy/2;

%             analytic(i,j) = -alpha*kappa/omega*E_0*cos(alpha*xi)*sin(beta*yj_stagger)*sin(kappa*omega*t);
            analytic(i,j) = E_0*sin(alpha*xi)*cos(beta*yj_stagger)*sin(kappa*omega*t);
        end
    end
%     analytic(1,:) = 0;
%     analytic(end,:) = 0;
end