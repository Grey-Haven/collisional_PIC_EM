function x = gauss_seidel(b,TOL)
    
    N = sqrt(length(b));

    x = zeros(length(b),1);
    MAX_ITERATIONS = 1e4;

    for iter = 1:MAX_ITERATIONS
        
        for i = 1:length(x)
            a_ii = -4;
            Ax_left = 0;
            Ax_rite = 0;
            Ax_lower = 0;
            Ax_upper = 0;
            if (i > 1 && mod(i-1,N) ~= 0)
                Ax_lower = x(i-1);
            end
            if (i < length(x) && mod(i,N) ~= 0)
                Ax_upper = x(i+1);
            end
            if (i - N >= 1)
                Ax_left = x(i-N);
            end
            if (i + N <= length(x))
                Ax_rite = x(i+N);
            end
            
            x(i) = 1/a_ii * (b(i) - Ax_left - Ax_rite - Ax_lower - Ax_upper);
        end
%         surf(reshape(x,N,N));
%         drawnow;
        
        Ax = laplacian_2D(x);
        residue = abs(norm(Ax - b));
        if (mod(iter,100) == 0)
            disp("STEP: " + iter + " RES: " + residue);
        end
        if (residue < TOL)
            break;
        end
    end
    
end