function [Ex_clean, Ey_clean] = divergence_clean_verb(Ex,Ey,rho,dx,dy,dt,sig_1)

%     Ex_clean = zeros(size(Ex));
%     Ey_clean = zeros(size(Ey));
    d = 1/(2*dt)*(dx*dy/(dx^2 + dy^2));

    Err = zeros(size(rho));
    Err(2:end-1,2:end-1) = dy*diff(Ey,1,2) + dx*diff(Ex,1,1) - rho(2:end-1,2:end-1)/sig_1;
    
    Ex_clean = Ex + dt/dx*(d*diff(Err(:,2:end-1),1,1));
    Ey_clean = Ey + dt/dy*(d*diff(Err(2:end-1,:),1,2));
    
end