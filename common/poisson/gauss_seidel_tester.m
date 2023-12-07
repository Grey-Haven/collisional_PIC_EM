clear;

N_iter = 5;

N = 4;

x_a = 0;
x_b = pi;
y_a = -pi/2;
y_b = pi/2;

dxs = zeros(N_iter,1);
errs = zeros(N_iter,1);

for iter = 1:N_iter

    dx = (x_b-x_a)/(N-1);
    dy = (y_b-y_a)/(N-1);

    x = x_a:dx:x_b;
    y = y_a:dy:y_b;

    RHS = -2*sin(x).*cos(y');    
    RHS_inner = RHS(2:end-1,2:end-1);
    b = reshape(RHS_inner,(N-2)^2,1);
    
    xi_approx_GS_vector = gauss_seidel(dx*dy*b,(dx*dy)^2);
    xi_approx_GS = reshape(xi_approx_GS_vector,N-2,N-2);

    xi_exact = sin(x(2:end-1)).*cos(y(2:end-1)');

    assert(max(max(abs(xi_exact-xi_approx_GS))) < dx*dy);
    
    errs(iter) = max(max(abs(xi_exact-xi_approx_GS)));
    dxs(iter) = dx;
    N = N*2;
end

linspace = dxs;
linear = linspace;
quadratic = linspace.^2;
cubic = linspace.^3;

loglog(dxs,errs);
% hold on;
% loglog(dxs,errs_GS);
% loglog(dxs,errs_sparse);
hold on;
loglog(linspace,linear);
loglog(linspace,quadratic);
loglog(linspace,cubic);
hold off;
legend("Errors", "GS Errors", "Sparse Errors", "linear", "quadratic", "cubic");
title("Refinement");
drawnow;


dx = (x_b-x_a)/(N-1);
dy = (y_b-y_a)/(N-1);

x = x_a:dx:x_b;
y = y_a:dy:y_b;

rho = zeros(N,N);
rho(floor(N/2),floor(N/2)) = 1;

rho_vector = reshape(rho(2:end-1,2:end-1),(N-2)^2,1);

TOL = dx*dy;
phi_GS_vec = dx*dy*gauss_seidel(rho_vector,TOL);
phi_approx = poisson_solver(rho(2:end-1,2:end-1),dx,dy);

A = construct_poisson_matrix(N-2);
% phi_GS_vec = reshape(phi_GS,(N-2)^2,1);
Aphi_GS_vec = A*phi_GS_vec;
assert(norm(Aphi_GS_vec - dx*dy*rho_vector) < (dx*dy)^2);