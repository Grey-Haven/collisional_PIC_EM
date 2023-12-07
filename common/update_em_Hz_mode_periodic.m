
tau = dt*10;
t_0 = 6*tau;

x_src_idx = Nx/4;
y_src_idx = Ny/4;

% periodic boundaries
Hz_down = Hz(:,end);
Hz_left = Hz(end,:);

Ex_up = Ex(:,1);
Ey_rite = Ey(1,:);
% =======================
% Update Curl of Hz field
% =======================
CEz = computeCurlE_z_rite_up_boundary(Ex,Ey,dx,dy,Ex_up,Ey_rite);

% =====================================
% Update H field from Curl of E fields
% =====================================
Hz = Hz + -dt.*CEz;
% pulse = energy_source(t, t_0, tau);
% Hz(x_src_idx,y_src_idx) = Hz(x_src_idx,y_src_idx) + pulse;

% ================================
% Compute Curl of Hx and Hy fields
% ================================
[CHx,CHy] = computeCurlH_xy_left_down_boundary(Hz,dx,dy,Hz_left,Hz_down);

% =====================================
% Update E fields from Curl of H field
% =====================================
Ex = Ex + dt*kappa^2*(CHx - sig_2*Jx);
Ey = Ey + dt*kappa^2*(CHy - sig_2*Jy);


function g_E = energy_source(t, t_0, tau)
    g_E = exp(-((t-t_0)/tau)^2);
end