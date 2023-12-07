N = 64;
A = construct_poisson_matrix(N);
x = rand(N*N,1);

Ax_analytic = A*x;
Ax_function = laplacian_2D(x);

Ax_a_r = reshape(Ax_analytic,N,N);
Ax_f_r = reshape(Ax_function,N,N);

disp(norm(Ax_analytic - Ax_function));