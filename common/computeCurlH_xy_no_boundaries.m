function [CHx,CHy] = computeCurlH_xy_no_boundaries(Hz, dx, dy)

    Nx = size(Hz,1);
    Ny = size(Hz,2);

    CHx = zeros(Nx,Ny-1);
    CHy = zeros(Nx-1,Ny);

    for nx = 1 : size(CHx,1)
        for ny = 1 : size(CHx,2)
            CHx(nx,ny) = (Hz(nx,ny+1) - Hz(nx,ny))/dy;
        end
    end
    for ny = 1 : size(CHy,2)
        for nx = 1 : size(CHy,1)
            CHy(nx,ny) = -(Hz(nx+1,ny) - Hz(nx,ny))/dx;
        end
    end
end