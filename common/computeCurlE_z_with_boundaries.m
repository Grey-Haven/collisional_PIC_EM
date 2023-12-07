function CEz = computeCurlE_z_with_boundaries(Ex, Ey, dx, dy, Ex_down, Ex_up, Ey_left, Ey_rite)

    [ExNx,ExNy] = size(Ex);
    [EyNx,EyNy] = size(Ey);
    Nx = max(ExNx, EyNx);
    Ny = max(ExNy, EyNy);
    CEz = zeros(Nx,Ny);
    
    Ex_with_boundaries = zeros(size(Ex,1),size(Ex,2)+2);
    Ex_with_boundaries(:,2:end-1) = Ex(:,:);
    Ex_with_boundaries(:,1) = Ex_down;
    Ex_with_boundaries(:,end) = Ex_up;
    
    Ey_with_boundaries = zeros(size(Ey,1)+2,size(Ey,2));
    Ey_with_boundaries(2:end-1,:) = Ey(:,:);
    Ey_with_boundaries(1,:) = Ey_left;
    Ey_with_boundaries(end,:) = Ey_rite;
    
    for ny = 1 : EyNy % larger of the two
        for nx = 1 : ExNx % larger of the two
            diffY = (Ey_with_boundaries(nx+1,ny) - Ey_with_boundaries(nx,ny))/dx;
            diffX = (Ex_with_boundaries(nx,ny+1) - Ex_with_boundaries(nx,ny))/dy;
            CEz(nx,ny) = diffY - diffX;
%             CEz(nx,ny) = (Ey_with_boundaries(nx+1,ny) - Ey_with_boundaries(nx,ny))/dx ...
%                        - (Ex_with_boundaries(nx,ny+1) - Ex_with_boundaries(nx,ny))/dy;
        end
    end
    
    % Ex_u/d should match the x index of Ex
    % Ey_l/r should match the y index of Ey
    % left/down should be subtracted
    % rite/up should be subtracting
        
%     CEz(1,1) = (Ey(1,1) - Ey_left(1))/dx - (Ex(1,1) - Ex_down(1))/dy;
%     for ny = 2 : ExNy % smaller of the two
%         CEz(1,ny) = (Ey(1,ny) - Ey_left(ny))/dx - (Ex(1,ny) - Ex(1,ny-1))/dy;
%         for nx = 2 : EyNx % smaller of the two
%             CEz(nx,ny) = (Ey(nx,ny) - Ey(nx-1,ny))/dx - (Ex(nx,ny) - Ex(nx,ny-1))/dy;
%         end
%         CEz(ExNx,ny) = (Ey_rite(ny) - Ey(ExNx-1,ny))/dx - (Ex(ExNx,ny) - Ex(ExNx,ny-1))/dy;
%     end
%     CEz(ExNx,1) = (Ey_rite(1) - Ey(ExNx-1,1))/dx - (Ex(ExNx,1) - Ex_down(ExNx))/dy;
% 
%     CEz(1,EyNy) = (Ey(1,EyNy) - Ey_left(EyNy))/dx - (Ex_up(1) - Ex(1,EyNy-1))/dy;    
%     for nx = 2 : EyNx
%         CEz(nx,1) = (Ey(nx,1) - Ey(nx-1,1))/dx - (Ex(nx,1) - Ex_down(nx))/dy;
%         CEz(nx,EyNy) = (Ey(nx,EyNy) - Ey(nx-1,EyNy))/dx - (Ex_up(nx) - Ex(nx,EyNy-1))/dy;
%     end
%     CEz(ExNx,EyNy) = (Ey_rite(EyNy) - Ey(ExNx-1,EyNy))/dx - (Ex_up(ExNx) - Ex(ExNx,EyNy-1))/dy;
%     
%     for i = 1:size(CEz,1)
%         for j = 1:size(CEz,2)
%             assert(CEz(i,j) == CEz_star(i,j));
%         end
%     end
    
end