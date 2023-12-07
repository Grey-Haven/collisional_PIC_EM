function [Jx,Jy] = scatter_particles_xy_periodic(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
                                        Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                                        x,y,dx,dy,dt)
    
    Jx = zeros(size(Jx,1), size(Jx,2));
    Jy = zeros(size(Jy,1), size(Jy,2));
    
    x_padded = x(1)-dx:dx:x(end)+dx;
    y_padded = y(1)-dy:dy:y(end)+dy;
    % We're periodic, so we need a place for the current at the edges to
    % land. We've padded x and y, however, the methods below will mod any
    % indices to place current on the proper node.
    Jx_padding_x = 1;
    Jx_padding_y = 1;
    Jy_padding_x = 1;
    Jy_padding_y = 1;
    cart_padding_x = 1;
    cart_padding_y = 1;
    scatter_particles_Jx_periodic;
    scatter_particles_Jy_periodic;
%     scatter_current_Jx_periodic;
%     scatter_current_Jy_periodic;
    
%     Jx_down = Jx(:,1);
%     Jx_up = Jx(:,end);
%     Jx_left = Jx(1,2:end-1);
%     Jx_rite = Jx(end,2:end-1);
%     
%     Jx(:,1) = Jx(:,1) + Jx_up;
%     Jx(:,end) = Jx(:,end) + Jx_down;
%     Jx(1,2:end-1) = Jx(1,2:end-1) + Jx_rite;
%     Jx(end,2:end-1) = Jx(end,2:end-1) + Jx_left;
%     
%     Jy_down = Jy(:,1);
%     Jy_up = Jy(:,end);
%     Jy_left = Jy(1,2:end-1);
%     Jy_rite = Jy(end,2:end-1);
%     
%     Jy(:,1) = Jy(:,1) + Jy_up;
%     Jy(:,end) = Jy(:,end) + Jy_down;
%     Jy(1,2:end-1) = Jy(1,2:end-1) + Jy_rite;
%     Jy(end,2:end-1) = Jy(end,2:end-1) + Jy_left;
    
    Jx = Jx ./ (dx*dy);
    Jy = Jy ./ (dx*dy);
%     Jx = Jx ./ (dx);
%     Jy = Jy ./ (dy);
end