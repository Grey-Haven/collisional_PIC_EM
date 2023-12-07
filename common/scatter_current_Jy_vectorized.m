
index_offset = 1;    
Jy_including_borders = zeros(size(Jy,1)+2,size(Jy,2));

locs = particles(:,1:3);
vels = particles(:,4:6);
tags = particles(:,7);
ms = particles(:,8);
ws = particles(:,9);
qs = particles(:,10);

mesh_padding = 1;

x_staggered = x;
y_staggered = y;

if (Jy_stagger_x ~= 0)
    x_staggered = x_padded(1:end-1)+Jy_stagger_x;
end

if (Jy_stagger_y ~= 0)
    y_staggered = y_padded(1:end-1)+Jy_stagger_y;
end

dt_orig = dt;x_prevs = prev_particle_locs(:,1);
y_prevs = prev_particle_locs(:,2);

x_currs = locs(:,1);
y_currs = locs(:,2);

% x_prevs_floored = max(0,x_prevs);
% x_currs_ceiled =  min(x(end),x_currs);
% Finding which node we are on (so no staggering)
lc_x_prevs = index_offset + (x_prevs+offset_x)/dx;
lc_y_prevs = index_offset + (y_prevs+offset_y)/dy;
lc_x_currs = index_offset + (x_currs+offset_x)/dx;
lc_y_currs = index_offset + (y_currs+offset_y)/dy;

i_currs = floor(lc_x_currs);
j_currs = floor(lc_y_currs);

i_prevs = floor(lc_x_prevs);
j_prevs = floor(lc_y_prevs);

% x_node_currs = x(i_currs);
% y_node_currs = y(j_currs);
% 
% x_node_prevs = x(i_prevs);
% y_node_prevs = y(j_prevs);

% particles_stay_within_node_idx = find(x_node_currs == x_node_prevs & y_node_currs == y_node_prevs);
% particles_leave_node_idx = find(~(x_node_currs == x_node_prevs & y_node_currs == y_node_prevs));

particles_stay_within_node_idx = find(i_currs == i_prevs & j_currs == j_prevs);
particles_leave_node_idx = find(~(i_currs == i_prevs & j_currs == j_prevs));

x_prevs_within_node = x_prevs(particles_stay_within_node_idx);
y_prevs_within_node = y_prevs(particles_stay_within_node_idx);

x_currs_within_node = x_currs(particles_stay_within_node_idx);
y_currs_within_node = y_currs(particles_stay_within_node_idx);

% x_node_currs_within_node = x_node_currs(particles_stay_within_node_idx);
% y_node_currs_within_node = y_node_currs(particles_stay_within_node_idx);
% 
% x_node_prevs_within_node = x_node_prevs(particles_stay_within_node_idx);
% y_node_prevs_within_node = y_node_prevs(particles_stay_within_node_idx);

qs_within_node = qs(particles_stay_within_node_idx);
ws_within_node = ws(particles_stay_within_node_idx);
        
[A1,A2] = scatter_I(x_prevs_within_node,y_prevs_within_node,x_currs_within_node,y_currs_within_node,...
                    offset_x,offset_y,...
                    0,0,0,0,...
                    qs_within_node,ws_within_node,x',y',dx,dy,dt);

Jy_lc_xs = Jy_padding_x+index_offset+(x_currs_within_node-Jy_stagger_x+offset_x)./dx;
Jy_lc_ys = Jy_padding_y+index_offset+(y_currs_within_node-Jy_stagger_y+offset_y)./dy;
Jy_is = floor(Jy_lc_xs);
Jy_js = round(Jy_lc_ys);

% idx1 = sub2ind(size(Jy_including_borders),Jy_is,Jy_js);
% idx2 = sub2ind(size(Jy_including_borders),Jy_is+1,Jy_js);

% Jy_including_borders(idx1) = Jy_including_borders(idx1) + A1;
% Jy_including_borders(idx2) = Jy_including_borders(idx2) + A2;
if (length(Jy_is) >= 1)
    Jy_including_borders_1 = accumarray([Jy_is,Jy_js],A1,size(Jy_including_borders));
    Jy_including_borders_2 = accumarray([Jy_is+1,Jy_js],A2,size(Jy_including_borders));
    Jy_including_borders = Jy_including_borders_1 + Jy_including_borders_2;
end
for p_i = 1:length(particles_leave_node_idx)
    
    idx = particles_leave_node_idx(p_i);
    
    q = qs(idx);
    w = ws(idx);
    x_curr = locs(idx,1);
    y_curr = locs(idx,2);
    
    x_prev = prev_particle_locs(idx,1);
    y_prev = prev_particle_locs(idx,2);
    
%     x_prev = max(0,x_prev); % questionable move here, logic is that there's no current outside the mesh anyways.
%     x_prev_floored = max(x(1),x_prev);
%     x_curr_ceiled =  min(x(end-1),x_curr);

    [x_prev_trim,y_prev_trim,x_curr_trim,y_curr_trim,del_t] = trim_particle_path(x_prev,y_prev,x_curr,y_curr,x,y,dt);

%     dt = dt_trim;
    x_prev = x_prev_trim;
    y_prev = y_prev_trim;
    x_curr = x_curr_trim;
    y_curr = y_curr_trim;

    % Finding which node we are on (so no staggering)
    lc_x_prev = index_offset + (x_prev_trim+offset_x)/dx;
    lc_y_prev = index_offset + (y_prev_trim+offset_y)/dy;
    lc_x_curr = index_offset + (x_curr_trim+offset_x)/dx;
    lc_y_curr = index_offset + (y_curr_trim+offset_y)/dy;
    
    i_curr = floor(lc_x_curr);
    j_curr = floor(lc_y_curr);
    
    i_prev = floor(lc_x_prev);
    j_prev = floor(lc_y_prev);
    
    if (x_prev == x(end))
        i_prev = i_prev - 1;
    end
    if (x_curr == x(end))
        i_curr = i_curr - 1;
    end
    
    if (y_prev == y(end))
        j_prev = j_prev - 1;
    end
    if (y_curr == y(end))
        j_curr = j_curr - 1;
    end
    
    x_node_curr = x(i_curr);
    y_node_curr = y(j_curr);
    
    x_node_prev = x(i_prev);
    y_node_prev = y(j_prev);
    
    if (x_node_curr == x_node_prev && y_node_curr == y_node_prev) % x and y are the same
        
        [line,line_inv] = line_from_points(x_prev,x_curr,y_prev,y_curr);
        
        prevIdxOffsetX = 0;
        prevIdxOffsetY = 0;
        currIdxOffsetX = 0;
        currIdxOffsetY = 0;

        if (x_curr == x(end))
            currIdxOffsetX = -1;
        end
        if (y_curr == y(end))
            currIdxOffsetY = -1;
        end

        [A3,A4] = scatter_I(x_prev,y_prev,x_curr,y_curr,offset_x,offset_y,...
                            prevIdxOffsetX,prevIdxOffsetY,currIdxOffsetX,currIdxOffsetY,...
                            q,w,x,y,dx,dy,dt);

        Jy_lc_x = Jy_padding_x+index_offset+((x_curr+x_prev)/2-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+((y_curr+y_prev)/2-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y);
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + A3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + A4;
        
    elseif (y_node_curr == y_node_prev) % y is the same
        
        m = (y_curr - y_prev)/(x_curr-x_prev);
        b = y_curr - m*x_curr;
        
        line = @(x_val) m*x_val+b;
        
        prevIdxOffsetX1 = 0;
        prevIdxOffsetY1 = 0;
        currIdxOffsetX1 = 0;
        currIdxOffsetY1 = 0;
        
        prevIdxOffsetX2 = 0;
        prevIdxOffsetY2 = 0;
        currIdxOffsetX2 = 0;
        currIdxOffsetY2 = 0;

        % Going from left to right
        if (x_node_curr > x_node_prev)
            x_mid = x(i_curr);
            currIdxOffsetX1 = -1;
        else % Going from right to left
            x_mid = x(i_prev);
            prevIdxOffsetX2 = -1;
        end

        if (x_curr == x(end))
            currIdxOffsetX2 = -1;
        end
        if (y_curr == y(end))
            currIdxOffsetY2 = -1;
        end

        y_mid = line(x_mid);
%         y_inter_end = line(x_inter_end);
        
        del_x = abs(x_curr - x_prev);
        del_y = abs(y_curr - y_prev);
        del_s = sqrt(del_x^2+del_y^2);
        
        del_x1 = abs(x_mid - x_prev);
        del_y1 = abs(y_mid - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        ft1 = (del_s1/del_s);
        ft2 = (1-ft1);
        
        del_t1 = dt*ft1;
        del_t2 = dt*ft2;
        
        assert(0 <= del_t1 && del_t1 <= 1);
        assert(0 <= del_t2 && del_t2 <= 1);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);

        [A3,A4] = scatter_I(x_prev,y_prev,x_mid,y_mid,offset_x,offset_y,...
                            prevIdxOffsetX1, prevIdxOffsetY1, currIdxOffsetX1, currIdxOffsetY1,...
                            q,w,x,y,dx,dy,dt);
        [B3,B4] = scatter_I(x_mid,y_mid,x_curr,y_curr,offset_x,offset_y,...
                            prevIdxOffsetX2, prevIdxOffsetY2, currIdxOffsetX2, currIdxOffsetY2,...
                            q,w,x,y,dx,dy,dt);
        
        lowerYIndex = 0;
        if (y_curr == y(end))
            lowerYIndex = -1;
        end

        Jy_lc_x = Jy_padding_x+index_offset+(x_prev-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+(y_prev-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y);
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + A3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + A4;

        Jy_lc_x = Jy_padding_x+index_offset+(x_curr-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+(y_curr-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y) + lowerYIndex;
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + B3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + B4;
        threshold = -1e-3;

    elseif (x_node_curr == x_node_prev) % x is the same
        
        [line, line_inv] = line_from_points(x_prev,x_curr,y_prev,y_curr);
        
        prevIdxOffsetX1 = 0;
        prevIdxOffsetY1 = 0;
        currIdxOffsetX1 = 0;
        currIdxOffsetY1 = 0;
        
        prevIdxOffsetX2 = 0;
        prevIdxOffsetY2 = 0;
        currIdxOffsetX2 = 0;
        currIdxOffsetY2 = 0;
        
        if (y_node_curr > y_node_prev)
            y_mid = y(j_curr);
            currIdxOffsetY1 = -1;
        else
            y_mid = y(j_prev);
            prevIdxOffsetY2 = -1;
        end

        if (x_curr == x(end))
            currIdxOffsetX2 = -1;
        end
        if (y_curr == y(end))
            currIdxOffsetY2 = -1;
        end
        
        x_mid = line_inv(y_mid);
        
        del_x = abs(x_curr - x_prev);
        del_y = abs(y_curr - y_prev);
        del_s = sqrt(del_x^2+del_y^2);
        
        del_x1 = abs(x_mid - x_prev);
        del_y1 = abs(y_mid - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        ft1 = (del_s1/del_s);
        ft2 = (1-ft1);
        
        del_t1 = dt*ft1;
        del_t2 = dt*ft2;
        
        assert(0 <= del_t1 && del_t1 <= 1);
        assert(0 <= del_t2 && del_t2 <= 1);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);

        [A3,A4] = scatter_I(x_prev,y_prev,x_mid,y_mid,offset_x,offset_y,...
                            prevIdxOffsetX1, prevIdxOffsetY1, currIdxOffsetX1, currIdxOffsetY1,...
                            q,w,x,y,dx,dy,dt);
        [B3,B4] = scatter_I(x_mid,y_mid,x_curr,y_curr,offset_x,offset_y,...
                            prevIdxOffsetX2, prevIdxOffsetY2, currIdxOffsetX2, currIdxOffsetY2,...
                            q,w,x,y,dx,dy,dt);
        lowerXIndex = 0;
        if (x_curr == x(end))
            lowerXIndex = -1;
        end
                        
        Jy_lc_x = Jy_padding_x+index_offset+(x_prev-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+(y_prev-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y);

        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + A3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + A4;

        Jy_lc_x = Jy_padding_x+index_offset+(x_curr-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+(y_curr-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x) + lowerXIndex;
        Jy_j = round(Jy_lc_y);
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + B3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + B4;
        
    else % Neither are the same
        
        [line, line_inv] = line_from_points(x_prev,x_curr,y_prev,y_curr);
        
        prevIdxOffsetX1 = 0;
        prevIdxOffsetY1 = 0;
        currIdxOffsetX1 = 0;
        currIdxOffsetY1 = 0;
        
        prevIdxOffsetX2 = 0;
        prevIdxOffsetY2 = 0;
        currIdxOffsetX2 = 0;
        currIdxOffsetY2 = 0;
        
        prevIdxOffsetX3 = 0;
        prevIdxOffsetY3 = 0;
        currIdxOffsetX3 = 0;
        currIdxOffsetY3 = 0;
        
        % left to right
        if (x_node_curr > x_node_prev)
            x_inter1 = x(i_curr);
            moveLeft = false;
        else % right to left
            x_inter1 = x(i_prev);
            moveLeft = true;
        end
        y_inter1 = line(x_inter1);
        
        % down to up
        if (y_node_curr > y_node_prev)
            y_inter2 = y(j_curr);
            moveDown = false;
        else % up to down
            y_inter2 = y(j_prev);
            moveDown = true;
        end
        x_inter2 = line_inv(y_inter2);
        
        % ie {x,y}_inter1 is closer to the starting point
        if (sqrt((x_inter1 - x_prev)^2 + (y_inter1 - y_prev)^2) < ...
            sqrt((x_inter2 - x_prev)^2 + (y_inter2 - y_prev)^2))
            x_mid1 = x_inter1;
            y_mid1 = y_inter1;
            
            x_mid2 = x_inter2;
            y_mid2 = y_inter2;

            if (moveLeft)
                prevIdxOffsetX2 = -1;
            else
                currIdxOffsetX1 = -1;
            end

            if (moveDown)
                prevIdxOffsetY3 = -1;
            else
                currIdxOffsetY2 = -1;
            end
        else
            x_mid1 = x_inter2; % Intentional number mixing (inter{1,2} doesn't indicate order)
            y_mid1 = y_inter2;
            
            x_mid2 = x_inter1;
            y_mid2 = y_inter1;

            if (moveDown)
                prevIdxOffsetY2 = -1;
            else
                currIdxOffsetY1 = -1;
            end
            
            if (moveLeft)
                prevIdxOffsetX3 = -1;
            else
                currIdxOffsetX2 = -1;
            end
        end

        if (x_curr == x(end))
            currIdxOffsetX3 = -1;
        end
        if (y_curr == y(end))
            currIdxOffsetY3 = -1;
        end
        
        del_x = abs(x_curr - x_prev);
        del_y = abs(y_curr - y_prev);
        del_s = sqrt(del_x^2+del_y^2);
        
        del_x1 = abs(x_mid1 - x_prev);
        del_y1 = abs(y_mid1 - y_prev);
        del_s1 = sqrt(del_x1^2+del_y1^2);
        
        del_x2 = abs(x_mid2 - x_mid1);
        del_y2 = abs(y_mid2 - y_mid1);
        del_s2 = sqrt(del_x2^2+del_y2^2);
        
        ft1 = (del_s1/del_s);
        ft2 = (del_s2/del_s);
        ft3 = (1-ft1-ft2);
        
        del_t1 = dt*ft1;
        del_t2 = dt*ft2;
        del_t3 = dt*ft3;
        
        assert(0 <= del_t1 && del_t1 <= 1);
        assert(0 <= del_t2 && del_t2 <= 1);
        assert(0 <= del_t3 && del_t3 <= 1);
        
        assert(0 <= ft1 && ft1 <= 1);
        assert(0 <= ft2 && ft2 <= 1);
        assert(0 <= ft3 && ft3 <= 1);
        
        [A3,A4] = scatter_I(x_prev,y_prev,x_mid1,y_mid1,offset_x,offset_y,...
                            prevIdxOffsetX1, prevIdxOffsetY1, currIdxOffsetX1, currIdxOffsetY1,...
                            q,w,x,y,dx,dy,dt);
        [B3,B4] = scatter_I(x_mid1,y_mid1,x_mid2,y_mid2,offset_x,offset_y,...
                            prevIdxOffsetX2, prevIdxOffsetY2, currIdxOffsetX2, currIdxOffsetY2,...
                            q,w,x,y,dx,dy,dt);
        [C3,C4] = scatter_I(x_mid2,y_mid2,x_curr,y_curr,offset_x,offset_y,...
                            prevIdxOffsetX3, prevIdxOffsetY3, currIdxOffsetX3, currIdxOffsetY3,...
                            q,w,x,y,dx,dy,dt);

        Jy_lc_x = Jy_padding_x+index_offset+(x_prev-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+(y_prev-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y);
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + A3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + A4;

        Jy_lc_x = Jy_padding_x+index_offset+((x_mid1+x_mid2)/2-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+((y_mid1+y_mid2)/2-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y);
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + B3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + B4;

        Jy_lc_x = Jy_padding_x+index_offset+(x_curr-Jy_stagger_x+offset_x)./dx;
        Jy_lc_y = Jy_padding_y+index_offset+(y_curr-Jy_stagger_y+offset_y)./dy;
        Jy_i = floor(Jy_lc_x);
        Jy_j = round(Jy_lc_y);
        
        Jy_including_borders(Jy_i,Jy_j) = Jy_including_borders(Jy_i,Jy_j) + C3;
        Jy_including_borders(Jy_i+1,Jy_j) = Jy_including_borders(Jy_i+1,Jy_j) + C4;
    end
    dt = dt_orig;
end

% for i = 1:size(Jy_including_borders,1)
%     for j = 1:size(Jy_including_borders,2)
%         if (Jy_including_borders(i,j)/dy < -2)
%             disp('foo');
%         end
%     end
% end

% x_nodes = x_staggered(is)';
% y_nodes = y_staggered(js)';
% 
% fxs = (locs(:,1) - x_nodes)./dx;
% fys = (locs(:,2) - y_nodes)./dy;

% Jy_including_borders(:,2) = Jy_including_borders(:,2) + Jy_including_borders(:,1);
% Jy_including_borders(:,end-1) = Jy_including_borders(:,end-1) + Jy_including_borders(:,end);
Jy(:,:) = Jy_including_borders(2:end-1,:); % cut out boundaries (perfectly conducting)

% Verboncoeur_2005_PIC_Review
function [Iy_1, Iy_2] = scatter_I(x_prev,y_prev,x_curr,y_curr,...
                                  offset_x,offset_y,...
                                  prevIdxOffsetX, prevIdxOffsetY, currIdxOffsetX, currIdxOffsetY,...
                                  q,w,x,y,dx,dy,dt)
    try
        index_offset = 1;

        lc_x_prev = (x_prev+offset_x)/dx; % Non staggered non padded grid
        lc_y_prev = (y_prev+offset_y)/dy;
        lc_x_curr = (x_curr+offset_x)/dx;
        lc_y_curr = (y_curr+offset_y)/dy;

        i_curr = floor(lc_x_curr) + index_offset + currIdxOffsetX;
        j_curr = floor(lc_y_curr) + index_offset + currIdxOffsetY;

        i_prev = floor(lc_x_prev) + index_offset + prevIdxOffsetX;
        j_prev = floor(lc_y_prev) + index_offset + prevIdxOffsetY;

        x_node_curr = x(i_curr);
        y_node_curr = y(j_curr);

        x_node_prev = x(i_prev);
        y_node_prev = y(j_prev);

    %     vx = vels(:,1);
    %     vy = vels(:,2);

        fx_curr = abs(x_curr - x_node_curr)/dx;
        fy_curr = abs(y_curr - y_node_curr)/dy;

        fx_prev = abs(x_prev - x_node_prev)/dx;
        fy_prev = abs(y_prev - y_node_prev)/dy;

        assert(all(fx_prev <= 1 & fx_prev >= 0));
        assert(all(fy_prev <= 1 & fy_prev >= 0));
        assert(all(fx_curr <= 1 & fx_curr >= 0));
        assert(all(fy_curr <= 1 & fy_curr >= 0));

        del_wy = fy_curr - fy_prev;
        wx_bar = (fx_curr + fx_prev)/2;

        Iy_1 = q./dt.*(1-wx_bar).*del_wy.*w;
        Iy_2 = q./dt.*wx_bar.*del_wy.*w;

%         Iy_1 = q*(1-wx_bar)*del_wy*w;
%         Iy_2 = q*wx_bar*del_wy*w;
    catch exception
        throw(exception)
    end
end