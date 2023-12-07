
tau = dt*10;
t_0 = 6*tau;

x_src_idx = Nx/2;
y_src_idx = Ny/2;

% Dirichlet boundaries of 0
Ex_down = zeros(length(Ex(:,1)),1);
Ex_up = zeros(length(Ex(:,end)),1);
Ey_left = zeros(length(Ey(1,:)),1);
Ey_rite = zeros(length(Ey(end,:)),1);

% ================================
% Compute Curl of Hx and Hy fields
% ================================
[CHx,CHy] = computeCurlH_xy_no_boundaries(Hz,dx,dy);

% =====================================
% Update E fields from Curl of H field
% =====================================
Ex = Ex + dt*kappa^2*(CHx - sig_2*Jx);
Ey = Ey + dt*kappa^2*(CHy - sig_2*Jy);

% =======================
% Update Curl of Hz field
% =======================
CEz = computeCurlE_z_with_boundaries(Ex,Ey,dx,dy,Ex_down,Ex_up,Ey_left,Ey_rite);

% =====================================
% Update H field from Curl of E fields
% =====================================
Hz = Hz + -dt.*CEz;
% pulse = energy_source(t, t_0, tau);
% Hz(x_src_idx,y_src_idx) = Hz(x_src_idx,y_src_idx) + pulse;


function g_E = energy_source(t, t_0, tau)
    g_E = exp(-((t-t_0)/tau)^2);
end