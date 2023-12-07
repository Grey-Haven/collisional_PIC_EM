function [Jx,Jy] = scatter_particles_xy(particles,prev_particle_locs,Jx,Jy,L_x,L_y,offset_x,offset_y,...
                                        Jx_stagger_x,Jx_stagger_y,Jy_stagger_x,Jy_stagger_y,...
                                        Jx_padding_x,Jx_padding_y,Jy_padding_x,Jy_padding_y,...
                                        x,y,dx,dy,dt)
    Jx = zeros(size(Jx,1)+Jx_padding_x*2, size(Jx,2)+Jx_padding_y*2);
    Jy = zeros(size(Jy,1)+Jy_padding_x*2, size(Jy,2)+Jy_padding_y*2);
    
    x_padded = x(1)-dx:dx:x(end)+dx;
    y_padded = y(1)-dy:dy:y(end)+dy;

%     [JxNx,JxNy] = size(Jx);
%     [JyNx,JyNy] = size(Jy);
%     scatter_particles_Jx;
%     scatter_particles_Jy;
    if (~isempty(particles))
        scatter_current_Jx_vectorized;
        scatter_current_Jy_vectorized;
%         scatter_particles_Jx;
%         scatter_particles_Jy;
    end
    
    Jx = Jx(1+Jx_padding_x:end-Jx_padding_x,1+Jx_padding_y:end-Jx_padding_y);
    Jy = Jy(1+Jy_padding_x:end-Jy_padding_x,1+Jy_padding_y:end-Jy_padding_y);
%     Jx = Jx ./ (dx*dy);
%     Jy = Jy ./ (dx*dy);
    Jx = Jx ./ (dx);
    Jy = Jy ./ (dy);
end