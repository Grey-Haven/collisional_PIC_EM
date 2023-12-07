function A = construct_poisson_matrix_sparse(N)
    K1D = spdiags(ones(N,1)*[1 -2 1],-1:1,N,N);  % 1d Poisson matrix
%     subplot(2,3,4), spy(K1D)

    I1D = speye(size(K1D));                       % 1d identity matrix
    A = kron(K1D,I1D)+kron(I1D,K1D);            % 2d PoissoN matrix
end