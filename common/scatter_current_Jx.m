
index_offset = 1;    
Jx_including_borders = zeros(size(Jx,1),size(Jx,2)+2);

locs = particles(:,1:3);
vels = particles(:,4:6);
tags = particles(:,7);
ms = particles(:,8);
ws = particles(:,9);
qs = particles(:,10);

% mesh_padding = 1;

x_staggered = x;
y_staggered = y;

if (Jx_stagger_x ~= 0)
    x_staggered = x_padded(1:end-1)+Jx_stagger_x;
end

if (Jx_stagger_y ~= 0)
    y_staggered = y_padded(1:end-1)+Jx_stagger_y;
end

eps = 1e-15;

dt_orig = dt;

for idx = 1:size(locs,1)
    q = qs(idx);
    w = ws(idx);
    x_curr = locs(idx,1);
    y_curr = locs(idx,2);
    
    x_prev = prev_particle_locs(idx,1);
    y_prev = prev_particle_locs(idx,2);

    x_prev_floored = max(0,x_prev);
    x_curr_ceiled =  min(x(end),x_curr)-eps;
    % Finding which node we are on (so no staggering)
    lc_x_prev = index_offset + (x_prev_floored+offset_x)/dx;
    lc_y_prev = index_offset + (y_prev+offset_y)/dy;
    lc_x_curr = index_offset + (x_curr_ceiled+offset_x)/dx;
    lc_y_curr = index_offset + (y_curr+offset_y)/dy;
    
    i_curr = floor(lc_x_curr);
    j_curr = floor(lc_y_curr);
    
    i_prev = floor(lc_x_prev);
    j_prev = floor(lc_y_prev);
    
    x_node_curr = x(i_curr);
    y_node_curr = y(j_curr);
    
    x_node_prev = x(i_prev);
    y_node_prev = y(j_prev);
    
    line = line_from_points(x_prev,x_curr,y_prev,y_curr);

    del_x = abs(x_curr - x_prev);
    del_y = abs(y_curr - y_prev);
    del_s = sqrt(del_x^2+del_y^2);

    del_x1 = abs(x_curr_ceiled - x_prev_floored);
    del_y1 = abs(y_curr - line(x_prev_floored));
    del_s1 = sqrt(del_x1^2+del_y1^2);

    ft1 = (del_s1/del_s);
    dt = ft1*dt_orig;

    x_prev = x_prev_floored;
    x_curr = x_curr_ceiled;
    
    if (x_node_curr == x_node_prev && y_node_curr == y_node_prev) % x and y are the same

        dt_sub = dt;
        x_prev_sub = x_prev;
        if (x_prev < 0)
            x_prev_sub = 0;
            m = (y_curr - y_prev)/(x_curr-x_prev);
            b = y_curr - m*x_curr;

            line = @(x_val) m*x_val+b;

            del_x = abs(x_curr - x_prev);
            del_y = abs(y_curr - y_prev);
            del_s = sqrt(del_x^2+del_y^2);

            del_x1 = abs(x_curr - x_prev_sub);
            del_y1 = abs(y_curr - line(x_prev_sub));
            del_s1 = sqrt(del_x1^2+del_y1^2);

            ft1 = (del_s1/del_s);
            dt_sub = ft1*dt;
        end
        
        [A1,A2] = scatter_I(x_prev_sub,y_prev,x_curr,y_curr,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,dt_sub);

        Jx_lc_x = Jx_padding_x+index_offset+(x_curr-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_curr-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + A1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + A2;

%         assert(A1 < 0);
%         assert(A2 < 0);
        
    elseif (y_node_curr == y_node_prev) % y is the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
        line = @(x_val) m*x_val+b;
        if (x_node_curr > x_node_prev)
            x_inter_end = x(i_curr)-eps;
            x_inter_beg = x(i_curr);
        else
            x_inter_end = x(i_prev);
            x_inter_beg = x(i_prev)-eps;
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
        
        assert(0 <= del_t1 && del_t1 <= 1);
        assert(0 <= del_t2 && del_t2 <= 1);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);
        
        [A1,A2] = scatter_I(x_prev,y_prev,x_inter_end,y_inter_end,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t1);
        [B1,B2] = scatter_I(x_inter_beg,y_inter_beg,x_curr,y_curr,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t2);

        Jx_lc_x = Jx_padding_x+index_offset+(x_prev-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_prev-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + A1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + A2;

        Jx_lc_x = Jx_padding_x+index_offset+(x_curr-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_curr-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + B1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + B2;
        
%         assert(A1 <= 0);
%         assert(A2 <= 0);
%         assert(B1 <= 0);
%         assert(B2 <= 0);
    elseif (x_node_curr == x_node_prev) % x is the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
%         line = @(x_val) m*x_val+b;
        line_inv = @(y_val) (y_val-b)/m;
        
        if (y_node_curr > y_node_prev)
            y_inter_beg = y(j_curr);
            y_inter_end = y(j_curr) - eps;
        else
            y_inter_beg = y(j_prev) - eps;
            y_inter_end = y(j_prev);
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
        
        assert(0 <= del_t1 && del_t1 <= 1);
        assert(0 <= del_t2 && del_t2 <= 1);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);

        [A1,A2] = scatter_I(x_prev,y_prev,x_inter_end,y_inter_end,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t1);
        [B1,B2] = scatter_I(x_inter_beg,y_inter_beg,x_curr,y_curr,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t2);

        Jx_lc_x = Jx_padding_x+index_offset+(x_prev-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_prev-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + A1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + A2;

        Jx_lc_x = Jx_padding_x+index_offset+(x_curr-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_curr-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + B1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + B2;
        
%         assert(A1 < 0);
%         assert(A2 < 0);
%         assert(B1 < 0);
%         assert(B2 < 0);

    else % Neither are the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
        line = @(x_val) m*x_val+b;
        line_inv = @(y_val) (y_val-b)/m;
        
        if (x_node_curr > x_node_prev)
            x_inter1_end = x(i_curr)-eps;
            x_inter1_beg = x(i_curr);
        else
            x_inter1_end = x(i_prev);
            x_inter1_beg = x(i_prev)-eps;
        end
        y_inter1_beg = line(x_inter1_beg);
        y_inter1_end = line(x_inter1_beg);
        
        if (y_node_curr > y_node_prev)
            y_inter2_beg = y(j_curr);
            y_inter2_end = y(j_curr) - eps;
        else
            y_inter2_beg = y(j_prev) - eps;
            y_inter2_end = y(j_prev);
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
        
        del_x1 = abs(x_inter_first_beg - x_prev);
        del_y1 = abs(y_inter_first_beg - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        del_x2 = abs(x_inter_second_beg - x_inter_first_beg);
        del_y2 = abs(y_inter_second_beg - y_inter_first_beg);
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
        
        assert(0 <= del_t1 && del_t1 <= 1);
        assert(0 <= del_t2 && del_t2 <= 1);
        assert(0 <= del_t3 && del_t3 <= 1);
        
        [A1,A2] = scatter_I(x_prev,y_prev,x_inter_first_end,y_inter_first_end,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t1);
        [B1,B2] = scatter_I(x_inter_first_beg,y_inter_first_beg,x_inter_second_end,y_inter_second_end,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t2);
        [C1,C2] = scatter_I(x_inter_second_beg,y_inter_second_beg,x_curr,y_curr,offset_x,offset_y,...
                            0,0,0,0,...
                            q,w,x,y,dx,dy,del_t3);

        Jx_lc_x = Jx_padding_x+index_offset+(x_prev-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_prev-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + A1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + A2;

        % Not sure if the epsilon makes a difference here.
        Jx_lc_x = Jx_padding_x+index_offset+(x_inter_first_beg-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_inter_first_beg-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + B1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + B2;

        Jx_lc_x = Jx_padding_x+index_offset+(x_curr-Jx_stagger_x+offset_x)./dx;
        Jx_lc_y = Jx_padding_y+index_offset+(y_curr-Jx_stagger_y+offset_y)./dy;
        Jx_i = round(Jx_lc_x);
        Jx_j = floor(Jx_lc_y);
        
        Jx_including_borders(Jx_i,Jx_j) = Jx_including_borders(Jx_i,Jx_j) + C1;
        Jx_including_borders(Jx_i,Jx_j+1) = Jx_including_borders(Jx_i,Jx_j+1) + C2;
%         assert(A1 < 0);
%         assert(A2 < 0);
%         assert(B1 < 0);
%         assert(B2 < 0);
%         assert(C1 < 0);
%         assert(C2 < 0);
    end
    dt = dt_orig;
end

% Jx_including_borders(2,:) = Jx_including_borders(2,:) + Jx_including_borders(1,:);
% Jx_including_borders(end-1,:) = Jx_including_borders(end-1,:) + Jx_including_borders(end,:);
Jx(:,:) = Jx_including_borders(:,2:end-1); % cut out boundaries (perfectly conducting)

% Verboncoeur_2005_PIC_Review
function [Ix_1, Ix_2] = scatter_I(x_prev,y_prev,x_curr,y_curr,...
                                  offset_x,offset_y,padding_x,padding_y,...
                                  stagger_x,stagger_y,q,w,x,y,dx,dy,dt)
    try
        index_offset = 1;

        lc_x_prev = (x_prev-stagger_x+offset_x)/dx; % Non staggered non padded grid
        lc_y_prev = (y_prev-stagger_y+offset_y)/dy;
        lc_x_curr = (x_curr-stagger_x+offset_x)/dx;
        lc_y_curr = (y_curr-stagger_y+offset_y)/dy;

        i_curr = floor(lc_x_curr) + padding_x + index_offset;
        j_curr = floor(lc_y_curr) + padding_y + index_offset;

        i_prev = floor(lc_x_prev) + padding_x + index_offset;
        j_prev = floor(lc_y_prev) + padding_y + index_offset;

        x_node_curr = x(i_curr);
        y_node_curr = y(j_curr);

        x_node_prev = x(i_prev);
        y_node_prev = y(j_prev);

    %     vx = vels(:,1);
    %     vy = vels(:,2);

        fx_curr = (x_curr - x_node_curr)/dx;
        fy_curr = (y_curr - y_node_curr)/dy;

        fx_prev = (x_prev - x_node_prev)/dx;
        fy_prev = (y_prev - y_node_prev)/dy;

        assert(fx_prev <= 1 && fx_prev >= 0);
        assert(fy_prev <= 1 && fy_prev >= 0);
        assert(fx_curr <= 1 && fx_curr >= 0);
        assert(fy_curr <= 1 && fy_curr >= 0);

        del_wx = fx_curr - fx_prev;
        wy_bar = (fy_curr + fy_prev)/2;

        Ix_1 = q/dt*del_wx*(1-wy_bar)*w;
        Ix_2 = q/dt*del_wx*wy_bar*w;

%         Ix_1 = q*del_wx*(1-wy_bar)*w;
%         Ix_2 = q*del_wx*wy_bar*w;
    catch exception
        throw(exception)
    end
end