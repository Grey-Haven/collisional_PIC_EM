disp(100*s/N_steps + "% complete.");

Ex_with_boundaries = zeros(size(Ex,1),size(Ex,2)+2);
Ex_with_boundaries(:,2:end-1) = Ex(:,:);
Ey_with_boundaries = zeros(size(Ey,1)+2,size(Ey,2));
Ey_with_boundaries(2:end-1,:) = Ey(:,:);

Jx_with_boundaries = zeros(size(Jx,1),size(Jx,2)+2);
Jx_with_boundaries(:,2:end-1) = Jx(:,:);
Jy_with_boundaries = zeros(size(Jy,1)+2,size(Jy,2));
Jy_with_boundaries(2:end-1,:) = Jy(:,:);

electrons = particles(particles(:,7) == eleTag,:);
ions = particles(particles(:,7) == ionTag,:);
ele_x = electrons(:,1);
ele_y = electrons(:,2);
ion_x = ions(:,1);
ion_y = ions(:,2);

vx_e = electrons(:,4);
vy_e = electrons(:,5);

var_vx = var(vx_e);
var_vy = var(vy_e);
var_v = (var_vx + var_vy)/2;
heat_hist((s/record_stride)+1) = var_v;
num_eles_hist((s/record_stride)+1) = length(electrons);

V0 = phi(floor(size(phi,1)), floor(size(phi,2)));
sheath_thickness = lambda_D*sqrt(2*V0/T_e);
s_hist((s/record_stride)+1) = sheath_thickness;
left = @(y) ones(length(y))*(x(1) + sheath_thickness);
rite = @(y) ones(length(y))*(x(end) - sheath_thickness);
bot = @(x) ones(length(x))*(y(1) + sheath_thickness);
top = @(x) ones(length(x))*(y(end) - sheath_thickness);

% histogram2(vx_e,vy_e,100);
% title("Velocity at t = " + t);
% 
% xlabel("v_x");
% ylabel("v_y");
% 
% xlim([-4,4]);
% ylim([-4,4]);
% zlim([0,350]);

        % scatter(ele_x,ele_y,'MarkerFaceColor',[0,1,1]);
        % hold on;
        % scatter(ion_x,ion_y,'MarkerFaceColor',[0,1,0],'MarkerEdgeAlpha',.1,'MarkerFaceAlpha',.1);
        % legend('Electrons','Ions');
        % hold off;
        % axis([a_x b_x a_y b_y]);
        % title("STEP: " + num2str(s) + ", t = " + num2str(t));
        % drawnow;

rho = scatter_rho(particles(:,1:3),particles(:,10),particles(:,9),x,y);
rho(1,:) = 0;
rho(end,:) = 0;
rho(:,1) = 0;
rho(:,end) = 0;

            subplot(2,2,1);
            surf(Ex_x,y,Ex_with_boundaries');
            shading interp;
            title("Ex");
            xlim([Ex_x(1),Ex_x(end)]);
            ylim([Ex_y(1),Ex_y(end)]);
            zlim([-1.5,1.5]);
            colorbar;
            clim([-1.5,1.5]);

            subplot(2,2,2);
            surf(x,Ey_y,Ey_with_boundaries');
            shading interp;
            title("Ey");
            xlim([Ey_x(1),Ey_x(end)]);
            ylim([Ey_y(1),Ey_y(end)]);
            zlim([-1.5,1.5]);
            colorbar;
            clim([-1.5,1.5]);

            subplot(2,2,3);
            phi_inner = -poisson_solver(rho(2:end-1,2:end-1)/sig_1,dx,dy);
            phi = zeros(size(phi_normal)+2);
            phi(2:end-1,2:end-1) = phi_inner;
            surf(x,y,phi);
            zlim([-1,8]);
            title("Phi");
    
            subplot(2,2,4);
            scatter(ele_x,ele_y,5,'MarkerFaceColor',[0,1,1]);
            hold on;
            axis([a_x b_x a_y b_y]);
            title('Electron Locations');
            hold off;
            sgtitle("Yee + Boris, " + col_title + ", "  + Nx + "x" + Ny + ", " + "STEP: " + num2str(s) + ", t = " + num2str(t));
            drawnow;

            set(gcf,'Position',[1250 250 1000 1000])
            hold off;

%         hold on;
%         scatter(particles(:,1), particles(:,2));
%         hold off;
currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);
%         
%         writeCSVfiles(Ex_with_boundaries',Ey_with_boundaries',...
%                       Jx_with_boundaries',Jy_with_boundaries',...
%                       Hz,...
%                       Ex_x,Ex_y,Ey_x,Ey_y,Hz_x,Hz_y,s,csvPath)