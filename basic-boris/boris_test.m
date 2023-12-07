% addpath(genpath([fileparts(pwd)]));
clear;
close all;

L_x = 20;
L_y = 20;
L_z = 20;

dx = .1;
dy = .1;
dz = .1;

kappa = 1;

Jx_stagger_x = 0;
Jx_stagger_y = 0;
Jy_stagger_x = 0;
Jy_stagger_y = 0;

Hz_stagger_x = 0;
Hz_stagger_y = 0;

Nx = L_x/dx + 1;
Ny = L_y/dy + 1;
Nz = L_z/dz + 1;

a_x = -10;
a_y = -10;
a_z = -10;
b_x = 10;
b_y = 10;
b_z = 10;

x = a_x:dx:b_x;
y = a_y:dy:b_y;
z = a_z:dz:b_z;

Hx_mag=0;
Hy_mag=0;
Hz_mag=0;

Ex_mag=0;
Ey_mag=0;
Ez_mag=0;

Hx = ones(Nx,Ny,Nz)*Hx_mag;
Hy = ones(Nx,Ny,Nz)*Hy_mag;
Hz = ones(Nx,Ny,Nz)*Hz_mag;

Ex = ones(Nx,Ny,Nz)*Ex_mag;
Ey = ones(Nx,Ny,Nz)*Ey_mag;
Ez = ones(Nx,Ny,Nz)*Ez_mag;

H = [Hx Hy Hz];
E = [Ex Ey Ez];

N_elec = 1;

eleTag = 45;
m_elec =  1;
q_elec = -1;
w = 1;
dt = 1/(kappa*sqrt(1/dx^2+1/dy^2));

particleX = [0,0,0];

ux = .967*kappa;
uy = 0;
uz = 0;

gammax = sqrt(1 + (ux/kappa)^2);
gammay = sqrt(1 + (uy/kappa)^2);
gammaz = sqrt(1 + (uz/kappa)^2);

particleU = [ux,uy,uz];

vx = ux;
vy = uy;
vz = uz;

particleV = [vx,vy,vz];
particles = [particleX(1),particleX(2),particleX(3),particleV(1),particleV(2),particleV(3),eleTag,m_elec,1,q_elec];

for i = 1:100
    boris_push;
    scatter3(particles(1),particles(2),particles(3),"red");
    hold on;
    drawnow;
end
% legend("Nonrelativistic");

vx = ux/gammax;
vy = uy/gammay;
vz = uz/gammaz;
particleV = [vx,vy,vz];

particles = [particleX(1),particleX(2),particleX(3),particleV(1),particleV(2),particleV(3),eleTag,m_elec,1,q_elec];
xs = [];
ys = [];
for i = 1:100
    boris_push_relativistic;
    xs = [xs,particles(1)];
    ys = [ys,particles(2)];
    scatter3(particles(1),particles(2),particles(3),"blue");
%     disp(norm(particles(4:6)));
    drawnow;
end
% legend("Relativistic");
xlabel("x");
ylabel("y");
zlabel("z");
title("Particle path in Hz field, initial velocity .967*kappa");

% B_mag = sqrt(Hx_mag^2+Hy_mag^2+Hz_mag^2);
% v_mag = norm(particleV);
% gamma = 1/(sqrt(1-(v_mag/kappa)^2));
% 
% omega_c = (q_elec*B_mag)/(gamma*m_elec); % cyclotron frequency
% dTheta = 2*atan(omega_c*dt/2);
% 
% omega = dTheta/dt;
% 
% x = [xs(1),xs(101),xs(201)];
% y = [ys(1),ys(101),ys(201)];
% 
% [c,r,x0,y0] = getCircle(x(1),y(1),x(2),y(2),x(3),y(3));
% 
% theta = linspace(0, 2*pi);
% x = r*cos(theta) + x0;
% y = r*sin(theta) + y0;
% plot(x, y);
% 
% R = sqrt(1 + (omega_c * dt/2)^2)*norm(particles(4:6)/omega_c);
% 
% diff = R-r;
% disp("R - r = " + diff);
% 
% diffs = zeros(length(x),1);
% 
% theta = acos(0/r - x0/r);
% for i = 1:350
%     x_i = r*cos(theta) + x0;
%     y_i = r*sin(theta) + y0;
%     theta = theta + dTheta;
%     scatter(x_i,y_i);
%     plot([x0,x_i],[y0,y_i]);
%     drawnow;
% end
% 
% for i = 3:length(x)
%     x_prev = x(i-2);
%     y_prev = y(i-2);
%     x_curr = x(i);
%     y_curr = y(i);
%     normDiff = sqrt((x_curr - x_prev)^2 + (y_curr - y_prev)^2);
%     diffs(i-2) = normDiff/2*sin(omega*dt) - r;
% %     disp(diffs(i-2));
% end

hold off;



function [c,r,x0,y0] = getCircle(x1,y1,x2,y2,x3,y3)    
    midpoint1 = [(x1 + x2)/2, (y1 + y2)/2];
    slope1 = (y2 - y1)/(x2 - x1);
    slope1_inv = -1/slope1;
    
    b1 = midpoint1(2) - midpoint1(1)*slope1;
    b1_inv = midpoint1(2) - midpoint1(1)*slope1_inv;
    
%     f1 = @(x) slope1*x+b1;
%     g1 = @(x) slope1_inv*x+b1_inv;
    
    % ie y = slope1_inv*x + b1
    
    midpoint2 = [(x2 + x3)/2, (y2 + y3)/2];
    slope2 = (y3 - y2)/(x3 - x2);
    slope2_inv = -1/slope2;
    
    b2 = midpoint2(2) - midpoint2(1)*slope2;
    b2_inv = midpoint2(2) - midpoint2(1)*slope2_inv;
    
%     f2 = @(x) slope2*x+b2;
%     g2 = @(x) slope2_inv*x+b2_inv;
    
    x0 = -1/(slope1_inv - slope2_inv) * (b1_inv - b2_inv);
    y0 = -1/(slope1_inv - slope2_inv) * (slope2_inv*b1_inv - slope1_inv*b2_inv);
    
%     min_x = min(x1, min(x2, x3));
%     max_x = max(x1, max(x2, x3));
    
%     x_star1 = min_x:.1:max_x;
%     x_star2 = min_x:.1:max_x;
    
%     scatter([x1,x2,x3],[y1,y2,y3]);
%     hold on;
%     scatter(midpoint1(1),midpoint1(2));
%     scatter(midpoint2(1),midpoint2(2));
%     plot(x_star1,f1(x_star1));
%     plot(x_star2,f2(x_star2));
%     plot(x_star1,g1(x_star1));
%     plot(x_star2,g2(x_star2));
%     scatter(x0,y0)
    
    % ie y = slope2_inv*x + b2
    c = [x0,y0];
    r = sqrt((x1 - x0)^2 + (y1 - y0)^2);
end