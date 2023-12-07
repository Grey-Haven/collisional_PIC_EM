function CEz = computeCurlE_z_rite_up_boundary(Ex, Ey, dx, dy, Ex_up, Ey_rite)

    [Nx,Ny] = size(Ex); % assumption, size(Ex) == size(Ey);
    CEz = zeros(Nx,Ny);
    
    Ex_with_boundaries = zeros(size(Ex,1),size(Ex,2)+1);
    Ex_with_boundaries(:,1:end-1) = Ex(:,:);
    Ex_with_boundaries(:,end) = Ex_up;
    
    Ey_with_boundaries = zeros(size(Ey,1)+1,size(Ey,2));
    Ey_with_boundaries(1:end-1,:) = Ey(:,:);
    Ey_with_boundaries(end,:) = Ey_rite;
    
    for ny = 1 : Ny
        for nx = 1 : Nx
            CEz(nx,ny) = (Ey_with_boundaries(nx+1,ny) - Ey_with_boundaries(nx,ny))/dx ...
                       - (Ex_with_boundaries(nx,ny+1) - Ex_with_boundaries(nx,ny))/dy;
        end
    end
    
%     for ny = 1 : Ny-1
%         for nx = 1 : Nx-1
%             CEz(nx,ny) = (Ey(nx+1,ny) - Ey(nx,ny))/dx ...
%                        - (Ex(nx,ny+1) - Ex(nx,ny))/dy;
%         end
%         CEz(Nx,ny) = (Ey( 1  ,ny) - Ey(Nx,ny))/dx ...
%                    - (Ex(Nx,ny+1) - Ex(Nx,ny))/dy;
%     end
%     for nx = 1:Nx-1
%         CEz(nx,Ny) = (Ey(nx+1,Ny) - Ey(nx,Ny))/dx ...
%                    - (Ex(nx, 1  ) - Ex(nx,Ny))/dy;
%     end
%     CEz(Nx,Ny) = (Ey(1,ny) - Ey(Nx,Ny))/dx ...
%                - (Ex(nx,1) - Ex(Nx,Ny))/dy;
    
end