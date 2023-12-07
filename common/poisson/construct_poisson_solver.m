function A_inv = construct_poisson_solver(N)
    A = construct_poisson_matrix(N);
    A_inv = inv(A);
end