function [CHx,CHy] = computeCurlH_xy_left_down_boundary(Hz, dx, dy, Hz_left, Hz_down)

    Nx = size(Hz,1);
    Ny = size(Hz,2);

    CHx = zeros(Nx,Ny); % assumption, size(Ex) == size(Ey);
    CHy = zeros(Nx,Ny);

%     for nx = 1 : Nx
%         CHx(nx,1) = (Hz(nx,1) - Hz_down(nx))/dy;
%         for ny = 2 : Ny
%             CHx(nx,ny) = (Hz(nx,ny) - Hz(nx,ny-1))/dy;
%         end
%     end
%     for ny = 1 : Ny
%         CHy(1,ny) = -(Hz(1,ny) - Hz_left(ny))/dx;
%         for nx = 2 : Nx
%             CHy(nx,ny) = -(Hz(nx,ny) - Hz(nx-1,ny))/dx;
%         end
%     end

    for nx = 1 : Nx
        CHx(nx,1) = (Hz(nx,1) - Hz(nx,end))/dy;
        for ny = 2 : Ny
            CHx(nx,ny) = (Hz(nx,ny) - Hz(nx,ny-1))/dy;
        end
    end
    for ny = 1 : Ny
        CHy(1,ny) = -(Hz(1,ny) - Hz(end,ny))/dx;
        for nx = 2 : Nx
            CHy(nx,ny) = -(Hz(nx,ny) - Hz(nx-1,ny))/dx;
        end
    end
end