function Ax = laplacian_2D(x)
    Ax = zeros(length(x),1);
    N = sqrt(length(x));
    for i = 1:length(x)
%         j = floor(i/N)+1;
%         if (mod(i-1,N) == 0 || mod(i,N) == 0 || j == 1 || j == N)
%             continue;
%         end
        Ax_left = 0;
        Ax_rite = 0;
        Ax_upper = 0;
        Ax_lower = 0;
        if (i > 1 && mod(i-1,N) ~= 0)
            Ax_lower = x(i-1);
        end
        if (i < length(x) && mod(i,N) ~= 0)
            Ax_upper = x(i+1);
        end
        if (i-N >= 1)
            Ax_left = x(i-N);
        end
        if (i+N <= length(x))
            Ax_rite = x(i+N);
        end
        Ax_center = -4*x(i);
        Ax(i) = Ax_left + Ax_rite + Ax_upper + Ax_lower + Ax_center;
    end
end