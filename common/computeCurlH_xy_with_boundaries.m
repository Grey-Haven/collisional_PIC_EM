function [CHx,CHy] = computeCurlH_xy_with_boundaries(Hz, dx, dy, Hz_down, Hz_up, Hz_left, Hz_rite)

    Nx = size(Hz,1);
    Ny = size(Hz,2);

    CHx = zeros(Nx,Ny+1);
    CHy = zeros(Nx+1,Ny);

    for nx = 1 : size(CHx,1)
        CHx(nx,1) = (Hz(nx,1) - Hz_down(nx))/dy;
        for ny = 2 : size(CHx,2)-1
            CHx(nx,ny) = (Hz(nx,ny) - Hz(nx,ny-1))/dy;
        end
        CHx(nx,Ny+1) = (Hz_up(nx) - Hz(nx,Ny))/dy;
    end
    for ny = 1 : size(CHy,2)
        CHy(1,ny) = -(Hz(1,ny) - Hz_left(ny))/dx;
        for nx = 2 : size(CHy,1)-1
            CHy(nx,ny) = -(Hz(nx,ny) - Hz(nx-1,ny))/dx;
        end
        CHy(Nx+1,ny) = -(Hz_rite(ny) - Hz(Nx,ny))/dx;
    end
end