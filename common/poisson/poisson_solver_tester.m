clear;

N_iter = 5;

N = 4;

x_a = 0;
x_b = pi;
y_a = -pi/2;
y_b = pi/2;

dxs = zeros(N_iter,1);
errs = zeros(N_iter,1);
% errs_GS = zeros(N_iter,1);
errs_sparse = zeros(N_iter,1);

for iter = 1:N_iter

    dx = (x_b-x_a)/(N-1);
    dy = (y_b-y_a)/(N-1);

    x = x_a:dx:x_b;
    y = y_a:dy:y_b;

    RHS = -2*sin(x).*cos(y');    
    RHS_inner = RHS(2:end-1,2:end-1);

    xi_approx = poisson_solver(RHS_inner,dx,dy);

    xi_exact = sin(x(2:end-1)).*cos(y(2:end-1)');

    assert(max(max(abs(xi_exact-xi_approx))) < dx*dy);
    
    errs(iter) = max(max(abs(xi_exact-xi_approx)));
    dxs(iter) = dx;
    N = N*2;
end

linspace = dxs;
linear = linspace;
quadratic = linspace.^2;
cubic = linspace.^3;

loglog(dxs,errs);
hold on;
loglog(linspace,linear);
loglog(linspace,quadratic);
loglog(linspace,cubic);
hold off;
legend("Errors", "linear", "quadratic", "cubic");
title("Refinement");
drawnow;