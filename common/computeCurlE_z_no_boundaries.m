function CEz = computeCurlE_z_no_boundaries(Ex, Ey, dx, dy)

    [ExNx,ExNy] = size(Ex);
    [EyNx,EyNy] = size(Ey);
    Nx = min(ExNx, EyNx);
    Ny = min(ExNy, EyNy);
    CEz = zeros(Nx,Ny);
    
    for ny = 1 : EyNy % smaller of the two
        for nx = 1 : ExNx % smaller of the two
            CEz(nx,ny) = (Ey(nx+1,ny) - Ey(nx,ny))/dx ...
                       - (Ex(nx,ny+1) - Ex(nx,ny))/dy;
        end
    end
end