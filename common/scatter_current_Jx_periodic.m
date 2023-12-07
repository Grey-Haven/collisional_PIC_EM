
index_offset = 1;

locs = particles(:,1:3);
vels = particles(:,4:6);
tags = particles(:,7);
ms = particles(:,8);
ws = particles(:,9);
qs = particles(:,10);

% mesh_padding = 1;

x_staggered = x_padded;
y_staggered = y_padded;

if (Jx_stagger_x ~= 0)
    x_staggered = x_padded(1:end-1)+Jx_stagger_x;
end

if (Jx_stagger_y ~= 0)
    y_staggered = y_padded(1:end-1)+Jx_stagger_y;
end

eps = 1e-15;

dt_orig = dt;

x_prevs = prev_particle_locs(:,1);
y_prevs = prev_particle_locs(:,2);

x_currs = locs(:,1);
y_currs = locs(:,2);

% Finding which node we are on (so no staggering)
lc_x_prevs = cart_padding_x + index_offset + (x_prevs+offset_x)/dx;
lc_y_prevs = cart_padding_y + index_offset + (y_prevs+offset_y)/dy;
lc_x_currs = cart_padding_x + index_offset + (x_currs+offset_x)/dx;
lc_y_currs = cart_padding_y + index_offset + (y_currs+offset_y)/dy;

i_currs = floor(lc_x_currs);
j_currs = floor(lc_y_currs);

i_prevs = floor(lc_x_prevs);
j_prevs = floor(lc_y_prevs);

x_node_currs = x_padded(i_currs);
y_node_currs = y_padded(j_currs);

x_node_prevs = x_padded(i_prevs);
y_node_prevs = y_padded(j_prevs);

particles_stay_within_node_idx = find(x_node_currs == x_node_prevs & y_node_currs == y_node_prevs);
particles_leave_node_idx = find(~(x_node_currs == x_node_prevs & y_node_currs == y_node_prevs));

x_prevs_within_node = x_prevs(particles_stay_within_node_idx);
y_prevs_within_node = y_prevs(particles_stay_within_node_idx);

x_currs_within_node = x_currs(particles_stay_within_node_idx);
y_currs_within_node = y_currs(particles_stay_within_node_idx);

x_node_currs_within_node = x_node_currs(particles_stay_within_node_idx);
y_node_currs_within_node = y_node_currs(particles_stay_within_node_idx);

x_node_prevs_within_node = x_node_prevs(particles_stay_within_node_idx);
y_node_prevs_within_node = y_node_prevs(particles_stay_within_node_idx);

qs_within_node = qs(particles_stay_within_node_idx);
ws_within_node = ws(particles_stay_within_node_idx);

[A1,A2] = scatter_I(x_prevs_within_node,y_prevs_within_node,x_currs_within_node,y_currs_within_node,offset_x,offset_y,...
                    0,0,0,0,...
                    qs_within_node,ws_within_node,x',y',dx,dy,dt);

Jx_lc_xs = physical_to_logical_space(x_currs_within_node,Jx_padding_x,Jx_stagger_x,offset_x,dx);
Jx_lc_ys = physical_to_logical_space(y_currs_within_node,Jx_padding_y,Jx_stagger_y,offset_y,dy);
Jx_is = round(Jx_lc_xs);
Jx_js = floor(Jx_lc_ys);

Jx_is = index_periodic(Jx_is,size(Jx,1)-1);
Jx_js = index_periodic(Jx_js,size(Jx,2)-1);

idx1 = sub2ind(size(Jx),Jx_is,Jx_js);
idx2 = sub2ind(size(Jx),Jx_is,index_periodic(Jx_js+1,size(Jx,2)-1));

Jx(idx1) = Jx(idx1) + A1;
Jx(idx2) = Jx(idx2) + A2;

for idx = particles_leave_node_idx
    
    dt = dt_orig;
    
    x_prev = x_prevs(idx);
    y_prev = y_prevs(idx);
    x_curr = x_currs(idx);
    y_curr = y_currs(idx);
    
    x_node_curr = x_node_currs(idx);
    x_node_prev = x_node_prevs(idx);
    y_node_curr = y_node_currs(idx);
    y_node_prev = y_node_prevs(idx);
    
    q = qs(idx);
    w = ws(idx);
    
    i_prev = i_prevs(idx);
    j_prev = j_prevs(idx);
    i_curr = i_currs(idx);
    j_curr = j_currs(idx);
    
    line = line_from_points(x_prev,x_curr,y_prev,y_curr);
    
    if (x_node_curr == x_node_prev && y_node_curr == y_node_prev) % x and y are the same
        
        [A1,A2] = scatter_I(x_prev_sub,y_prev,x_curr,y_curr,offset_x,offset_y,...
                            cart_padding_x,cart_padding_y,0,0,...
                            q,w,x_padded,y_padded,dx,dy,dt);

        Jx_lc_x = physical_to_logical_space(x_curr,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_curr,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + A1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + A2;
        
%         assert(~isnan(A1));
%         assert(~isnan(A2));
        
    elseif (y_node_curr == y_node_prev) % y is the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
        line = @(x_val) m*x_val+b;
        if (x_node_curr > x_node_prev)
            x_inter_end = x_padded(i_curr)-eps;
            x_inter_beg = x_padded(i_curr);
        else
            x_inter_end = x_padded(i_prev);
            x_inter_beg = x_padded(i_prev)-eps;
        end
        y_inter_beg = line(x_inter_beg);
        y_inter_end = line(x_inter_end);
        
        del_x = abs(x_curr - x_prev);
        del_y = abs(y_curr - y_prev);
        del_s = sqrt(del_x^2+del_y^2);
        
        del_x1 = abs(x_inter_end - x_prev);
        del_y1 = abs(y_inter_end - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        ft1 = (del_s1/del_s);
        ft2 = (1-ft1);
        
        del_t1 = dt*ft1;
        del_t2 = dt*ft2;
        
        assert(0 <= del_t1 && del_t1 <= dt);
        assert(0 <= del_t2 && del_t2 <= dt);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);

        if (del_t1 == 0)
            A1 = 0;
            A2 = 0;
        else
            [A1,A2] = scatter_I(x_prev,y_prev,x_inter_end,y_inter_end,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end
        if (del_t2 == 0)
            B1 = 0;
            B2 = 0;
        else
            [B1,B2] = scatter_I(x_inter_beg,y_inter_beg,x_curr,y_curr,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end

        Jx_lc_x = physical_to_logical_space(x_prev,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_prev,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + A1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + A2;

        Jx_lc_x = physical_to_logical_space(x_curr,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_curr,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + B1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + B2;
        
%         assert(~isnan(A1));
%         assert(~isnan(A2));
%         assert(~isnan(B1));
%         assert(~isnan(B2));

    elseif (x_node_curr == x_node_prev) % x is the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
%         line = @(x_val) m*x_val+b;
        line_inv = @(y_val) (y_val-b)/m;
        
        if (y_node_curr > y_node_prev)
            y_inter_beg = y_padded(j_curr);
            y_inter_end = y_padded(j_curr) - eps;
        else
            y_inter_beg = y_padded(j_prev) - eps;
            y_inter_end = y_padded(j_prev);
        end
        x_inter_beg = line_inv(y_inter_beg);
        x_inter_end = line_inv(y_inter_end);
        
        del_x = abs(x_curr - x_prev);
        del_y = abs(y_curr - y_prev);
        del_s = sqrt(del_x^2+del_y^2);
        
        del_x1 = abs(x_inter_end - x_prev);
        del_y1 = abs(y_inter_end - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        ft1 = (del_s1/del_s);
        ft2 = (1-ft1);
        
        del_t1 = dt*ft1;
        del_t2 = dt*ft2;
        
        assert(0 <= del_t1 && del_t1 <= dt);
        assert(0 <= del_t2 && del_t2 <= dt);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);

        if (del_t1 == 0)
            A1 = 0;
            A2 = 0;
        else
            [A1,A2] = scatter_I(x_prev,y_prev,x_inter_end,y_inter_end,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end
        if (del_t2 == 0)
            B1 = 0;
            B2 = 0;
        else
            [B1,B2] = scatter_I(x_inter_beg,y_inter_beg,x_curr,y_curr,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end

        Jx_lc_x = physical_to_logical_space(x_prev,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_prev,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + A1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + A2;

        Jx_lc_x = physical_to_logical_space(x_curr,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_curr,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + B1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + B2;
        
%         assert(~isnan(A1));
%         assert(~isnan(A2));
%         assert(~isnan(B1));
%         assert(~isnan(B2));

    else % Neither are the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
        line = @(x_val) m*x_val+b;
        line_inv = @(y_val) (y_val-b)/m;
        
        if (x_node_curr > x_node_prev)
            x_inter1_end = x_padded(i_curr)-eps;
            x_inter1_beg = x_padded(i_curr);
        else
            x_inter1_end = x_padded(i_prev);
            x_inter1_beg = x_padded(i_prev)-eps;
        end
        y_inter1_beg = line(x_inter1_beg);
        y_inter1_end = line(x_inter1_beg);
        
        if (y_node_curr > y_node_prev)
            y_inter2_beg = y_padded(j_curr);
            y_inter2_end = y_padded(j_curr) - eps;
        else
            y_inter2_beg = y_padded(j_prev) - eps;
            y_inter2_end = y_padded(j_prev);
        end
        x_inter2_beg = line_inv(y_inter2_beg);
        x_inter2_end = line_inv(y_inter2_end);
        
        if (sqrt((x_inter1_beg-x_prev)^2 + (y_inter1_beg-y_prev)^2) < ...
            sqrt((x_inter2_beg-x_prev)^2 + (y_inter2_beg-y_prev)^2))
            x_inter_first_beg = x_inter1_beg;
            x_inter_first_end = x_inter1_end;
            y_inter_first_beg = y_inter1_beg;
            y_inter_first_end = y_inter1_end;
            
            x_inter_second_beg = x_inter2_beg;
            x_inter_second_end = x_inter2_end;
            y_inter_second_beg = y_inter2_beg;
            y_inter_second_end = y_inter2_end;
        else
            x_inter_first_beg = x_inter2_beg;
            x_inter_first_end = x_inter2_end;
            y_inter_first_beg = y_inter2_beg;
            y_inter_first_end = y_inter2_end;
            
            x_inter_second_beg = x_inter1_beg;
            x_inter_second_end = x_inter1_end;
            y_inter_second_beg = y_inter1_beg;
            y_inter_second_end = y_inter1_end;
        end
        
        del_x = abs(x_curr - x_prev);
        del_y = abs(y_curr - y_prev);
        del_s = sqrt(del_x^2+del_y^2);
        
        del_x1 = abs(x_inter_first_end - x_prev);
        del_y1 = abs(y_inter_first_end - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        del_x2 = abs(x_inter_second_end - x_inter_second_beg);
        del_y2 = abs(y_inter_second_end - y_inter_second_beg);
        del_s2 = sqrt(del_x2^2+del_y2^2);
        
        ft1 = (del_s1/del_s);
        ft2 = (del_s2/del_s);
        ft3 = (1-ft1-ft2);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);
        assert(0 <= ft3 && ft3 <= 1);
        
        del_t1 = dt*ft1;
        del_t2 = dt*ft2;
        del_t3 = dt*ft3;
        
        assert(0 <= del_t1 && del_t1 <= dt);
        assert(0 <= del_t2 && del_t2 <= dt);
        assert(0 <= del_t3 && del_t3 <= dt);

        if (del_t1 == 0)
            A1 = 0;
            A2 = 0;
        else
            [A1,A2] = scatter_I(x_prev,y_prev,x_inter_first_end,y_inter_first_end,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end
        if (del_t2 == 0)
            B1 = 0;
            B2 = 0;
        else
            [B1,B2] = scatter_I(x_inter_first_beg,y_inter_first_beg,x_inter_second_end,y_inter_second_end,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end
        if (del_t3 == 0)
            C1 = 0;
            C2 = 0;
        else
            [C1,C2] = scatter_I(x_inter_second_beg,y_inter_second_beg,x_curr,y_curr,offset_x,offset_y,...
                                cart_padding_x,cart_padding_y,0,0,...
                                q,w,x_padded,y_padded,dx,dy,dt);
        end

        Jx_lc_x = physical_to_logical_space(x_prev,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_prev,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + A1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + A2;

        Jx_lc_x = physical_to_logical_space(x_inter_first_beg,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_inter_first_beg,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + B1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + B2;

        Jx_lc_x = physical_to_logical_space(x_curr,Jx_padding_x,Jx_stagger_x,offset_x,dx);
        Jx_lc_y = physical_to_logical_space(y_curr,Jx_padding_y,Jx_stagger_y,offset_y,dy);
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_i = index_periodic(Jx_i,size(Jx,1)-1);
        Jx_j = index_periodic(Jx_j,size(Jx,2)-1);
        
        Jx(Jx_i,Jx_j) = Jx(Jx_i,Jx_j) + C1;
        Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) = Jx(Jx_i,index_periodic(Jx_j+1,size(Jx,2)-1)) + C2;
        
%         assert(~isnan(A1));
%         assert(~isnan(A2));
%         assert(~isnan(B1));
%         assert(~isnan(B2));
%         assert(~isnan(C1));
%         assert(~isnan(C2));

    end
end



% Verboncoeur_2005_PIC_Review
function [Ix_1, Ix_2] = scatter_I(x_prev,y_prev,x_curr,y_curr,...
                                  offset_x,offset_y,padding_x,padding_y,...
                                  stagger_x,stagger_y,q,w,x,y,dx,dy,dt)
    try
        index_offset = 1;

        lc_x_prev = padding_x + index_offset + (x_prev-stagger_x+offset_x)/dx; % Non staggered non padded grid
        lc_y_prev = padding_y + index_offset + (y_prev-stagger_y+offset_y)/dy;
        lc_x_curr = padding_x + index_offset + (x_curr-stagger_x+offset_x)/dx;
        lc_y_curr = padding_y + index_offset + (y_curr-stagger_y+offset_y)/dy;

        i_curr = floor(lc_x_curr);
        j_curr = floor(lc_y_curr);

        i_prev = floor(lc_x_prev);
        j_prev = floor(lc_y_prev);

        x_node_curr = x(i_curr);
        y_node_curr = y(j_curr);

        x_node_prev = x(i_prev);
        y_node_prev = y(j_prev);

    %     vx = vels(:,1);
    %     vy = vels(:,2);

        wx_curr = (x_curr - x_node_curr);
        wy_curr = (y_curr - y_node_curr);

        wx_prev = (x_prev - x_node_prev);
        wy_prev = (y_prev - y_node_prev);

        assert(all(wx_prev <= dx & wx_prev >= 0));
        assert(all(wy_prev <= dy & wy_prev >= 0));
        assert(all(wx_curr <= dx & wx_curr >= 0));
        assert(all(wy_curr <= dy & wy_curr >= 0));

        del_wx = wx_curr - wx_prev;
        wy_bar = (wy_curr + wy_prev)/2;

        Ix_1 = q./dt.*del_wx.*(1-wy_bar).*w;
        Ix_2 = q./dt.*del_wx.*wy_bar.*w;

%         Ix_1 = q*del_wx*(1-wy_bar)*w;
%         Ix_2 = q*del_wx*wy_bar*w;
    catch exception
        throw(exception)
    end
end