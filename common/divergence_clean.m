function [Ex_clean, Ey_clean] = divergence_clean(Ex,Ey,rho,dx,dy,sig_1,A_inv)

    Ex_x = diff(Ex,1,1)/dx;
    Ey_y = diff(Ey,1,2)/dy;
    
    div_E = Ex_x + Ey_y;

    RHS = div_E - rho(2:end-1,2:end-1)/sig_1;
    xi = zeros(size(rho,1),size(rho,2));
    if ~(exist('A_inv','var'))
        xi(2:end-1,2:end-1) = poisson_solver(RHS,dx,dy);
    else
        xi(2:end-1,2:end-1) = poisson_solver(RHS,dx,dy,A_inv);
    end
    xi_x = diff(xi,1,1)/dx;
    xi_y = diff(xi,1,2)/dx;

    Ex_clean = Ex - xi_x(:,2:end-1);
    Ey_clean = Ey - xi_y(2:end-1,:);

%     xi_fig = figure;
%     surf(xi);
%     title("xi");
%     
%     xi_x_fig = figure;
%     surf(xi_x);
%     title("xi_x");
%     
%     Ex_fig = figure;
%     surf(Ex);
%     title("Ex");
%     
%     Ex_x_fig = figure;
%     surf(Ex_x);
%     title("Ex_x");
    
%     res = Ex_x(:,2:end-1) + Ey_y(2:end-1,:) - rho(2:end-1,2:end-1)/sig_1;
%     surf(res);
%     drawnow;

end