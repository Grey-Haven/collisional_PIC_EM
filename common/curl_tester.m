% [CHx, CHy] = computeCurlH_xy(Hz, dx, dy);
clear;
close all;
iters = 5;

dxs = zeros(1,iters);
diffs = zeros(1,iters);

for d = 0:iters-1

    dx = .1 * 2^-d;
    dy = .1 * 2^-d;

    x = 0:dx:1-dx;
    y = 0:dy:1-dy;
    
    Ex_x = x + dx/2;
    Ex_y = y;
    
    Ey_x = x;
    Ey_y = y + dy/2;

    Ex = Ex_function(Ex_x',Ex_y);
    Ey = Ey_function(Ey_x',Ey_y);

    Ex_up = Ex(:,1);
    Ey_rite = Ey(1,:);

    CEz_numeric = computeCurlE_z_rite_up_boundary(Ex, Ey, dx, dy, Ex_up, Ey_rite);

    CEz_analytic = Ey_dx_function(Ex_x',Ey_y) - Ex_dy_function(Ex_x',Ey_y);

    CEz_analytic_loop = zeros(size(CEz_analytic));

    for i = 1:size(CEz_analytic_loop,1)
        for j = 1:size(CEz_analytic_loop,2)
            x_i = x(i);
            y_j = y(j);
            CEz_analytic_loop(i,j) = Ey_dx_function(x_i,y_j) ...
                                   - Ex_dy_function(x_i,y_j);
        end
    end

    diff = abs(CEz_numeric - CEz_analytic);

    subplot(4,1,1);
    surf(Ex_x,Ey_y,CEz_numeric);
    title("Numerical");

    subplot(4,1,2);
    surf(Ex_x,Ey_y,CEz_analytic);
    title("Analytic");

    subplot(4,1,3);
    surf(Ex_x,Ey_y,CEz_analytic_loop);
    title("Loop");

    subplot(4,1,4);
    surf(Ex_x,Ey_y,diff);
    title("Abs Diff");
    drawnow;
    
    dxs(d+1) = dx;
    diffs(d+1) = max(max(abs(CEz_numeric - CEz_analytic)));
end

figure;
loglog(dxs,diffs);
linspace = dxs;
linear = linspace;
quadratic = linspace.^2;
cubic = linspace.^3;
hold on;
loglog(linspace,linear);
loglog(linspace,quadratic);
loglog(linspace,cubic);
legend("Errors", "linear", "quadratic", "cubic");


function Ex = Ex_function(x,y)
    Ex = sin(2*pi*x).*sin(2*pi*y);
end

function Ex_dy = Ex_dy_function(x,y)
    Ex_dy = 2*pi*sin(2*pi*x).*cos(2*pi*y);
end

function Ey = Ey_function(x,y)
    Ey = sin(4*pi*x).*sin(4*pi*y);
end

function Ey_dx = Ey_dx_function(x,y)
    Ey_dx = 4*pi*cos(4*pi*x).*sin(4*pi*y);
end

% dirichlet
% function Ex = Ex_function(x,y)
%     Ex = x.*y.*(x-1).*(y-1);
% end
% 
% function Ex_dy = Ex_dy_function(x,y)
%     Ex_dy = x.*(x-1).*(2*y-1);
% end
% 
% function Ey = Ey_function(x,y)
%     Ey = x.^2.*y.^2.*(x-1).*(y-1);
% end
% 
% function Ey_dx = Ey_dx_function(x,y)
%     Ey_dx = x.*y.^2.*(3*x-2).*(y-1);
% end