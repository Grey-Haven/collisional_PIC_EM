function A = construct_poisson_matrix(N)
    A = zeros(N*N,N*N);
    T = diag(-4*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1);
    I = -eye(N,N);
    A(1:N,1:N) = T;
    A(1:N,N+1:N+N) = I;
    for i = 2:N-1
        A((i-1)*N+1:(i)*N,(i-2)*N+1:(i-1)*N) = I;
        A((i-1)*N+1:(i)*N,(i-1)*N+1:(i)*N) = T;
        A((i-1)*N+1:(i)*N,(i)*N+1:(i+1)*N) = I;
    end
    A(N*N-N+1:N*N,N*N-N-N+1:N*N-N) = I;
    A(N*N-N+1:N*N,N*N-N+1:N*N) = T;
end