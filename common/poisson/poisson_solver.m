function xi = poisson_solver(RHS,dx,dy,preconstructed_solver)
    assert(size(RHS,1)==size(RHS,2)); % designed for a perfectly square box
    N = size(RHS,1);
    % Vectorize the rhs
    RHS_vector = zeros(N*N,1);
    for row = 1:N
        for col = 1:N
            RHS_vector(N*(row-1) + col) = RHS(row,col);
        end
    end
    % Construct the 5 point FD matrix A for the system A*xi = RHS
    % -1,...,-1,4,-1,...,-1
    if ~(exist('preconstructed_solver','var'))
        A = construct_poisson_matrix_sparse(N);
%         A_inv = construct_poisson_solver(N);
%         xi_vector = A_inv*(dx*dy*RHS_vector);
        xi_vector = A\(dx*dy*RHS_vector);
    else
        xi_vector = preconstructed_solver*(dx*dy*RHS_vector);
    end
    
    % un-vectorize the solution
    xi = zeros(N,N);
    for row = 1:N
        for col = 1:N
            xi(row,col) = xi_vector(N*(row-1)+ col);
        end
    end
end