% try
    index_offset = 1;
    Jy_including_borders = zeros(size(Jy,1)+2,size(Jy,2));

    x_staggered = x;
    y_staggered = y;

    if (Jy_stagger_x ~= 0)
        x_staggered = x_padded(1:end-1)+Jy_stagger_x;
    end

    if (Jy_stagger_y ~= 0)
        y_staggered = y_padded(1:end-1)+Jy_stagger_y;
    end

    locs = particles(:,1:3);
    vels = particles(:,4:6);
    tags = particles(:,7);
    ms = particles(:,8);
    ws = particles(:,9);
    qs = particles(:,10);
%     ws = weightMap(char(tags));

    lc_x = Jy_padding_x+index_offset+(locs(:,1)-Jy_stagger_x+offset_x)./dx;
    lc_y = Jy_padding_y+index_offset+(locs(:,2)-Jy_stagger_y+offset_y)./dy;
    is = floor(lc_x);
    js = floor(lc_y);
    wys = qs.*ws.*vels(:,2);
    
    x_nodes = x_staggered(is)';
    y_nodes = y_staggered(js)';

    fxs = (locs(:,1) - x_nodes)./dx;
    fys = (locs(:,2) - y_nodes)./dy;
    
    Jy5 = accumarray([is,js],(1-fxs).*(1-fys).*wys,size(Jy_including_borders));
    Jy6 = accumarray([is,js+1],(1-fxs).*fys.*wys,size(Jy_including_borders));
    Jy7 = accumarray([is+1,js],fxs.*(1-fys).*wys,size(Jy_including_borders));
    Jy8 = accumarray([is+1,js+1],fxs.*fys.*wys,size(Jy_including_borders));

    Jy_including_borders_vec2 = Jy5+Jy6+Jy7+Jy8;
    Jy_including_borders = Jy_including_borders_vec2;
    

%     for p = 1:size(particles,1)
%         particle = particles(p,:);
%         
%         locs = particle(1:3);
%         vels = particle(4:6);
% %         tags = particle(7);
%         weight = particle(9);
%         q = particle(10);
%         
%         lc_x = Jy_padding_x+index_offset+(locs(1)-Jy_stagger_x+offset_x)/dx;
%         lc_y = Jy_padding_y+index_offset+(locs(2)-Jy_stagger_y+offset_y)/dy;
%         
%         i = floor(lc_x);
%         j = floor(lc_y);
%         
% %         q = particle.charge;
%         
%         vy = vels(2);
% %         weight = weightMap(particle.tag);
%         wy = q*vy*weight;
% 
% %         x_staggered = x;
% %         y_staggered = y;
% 
% %         if (Jy_stagger_x ~= 0)
% %             x_staggered = x_padded(1:end-1)+Jy_stagger_x;
% %         end
% 
% %         if (Jy_stagger_y ~= 0)
% %         y_staggered = y_padded(1:end-1)+Jy_stagger_y;
% %         end
%         
%         x_node = x_staggered(i); % - Jy_stagger_x;
%         y_node = y_staggered(j); % - Jy_stagger_y;
% 
%         fx = (locs(1) - x_node)/dx;
%         fy = (locs(2) - y_node)/dy;
%         
%         if (fx < 0 || fy < 0 || fx > 1 || fy > 1)
%             throw(exception);
%         end
% 
%         Jy_including_borders(i,j) = Jy_including_borders(i,j) + (1-fx)*(1-fy)*wy;
%         Jy_including_borders(i,j+1) = Jy_including_borders(i,j+1) + (1-fx)*fy*wy;
%         Jy_including_borders(i+1,j) = Jy_including_borders(i+1,j) + fx*(1-fy)*wy;
%         Jy_including_borders(i+1,j+1) = Jy_including_borders(i+1,j+1) + fx*fy*wy;
%         
%     end
%     if ~(abs(norm(Jy_including_borders-Jy_including_borders_vec2)) < eps)
%         disp('foo');
%     end
%     assert(isequal(Jy_including_borders,Jy_including_borders_vec2));
    Jy_including_borders(:,2) = Jy_including_borders(:,2) + Jy_including_borders(:,1);
    Jy_including_borders(:,end-1) = Jy_including_borders(:,end-1) + Jy_including_borders(:,end);
    Jy(:,:) = Jy_including_borders(2:end-1,:);
%     Jy(:,2) = Jy(:,2) + Jy(:,1);
%     Jy(:,end-1) = Jy(:,end-1) + Jy(:,end);
%     if (Jy_padding_x == 1)
%         % TRIAL
%         Jy(2,:) = Jy(2,:) + Jy(1,:);
%         Jy(end-1,:) = Jy(end-1,:) + Jy(end,:);
%         % END TRIAL
%         Jy(1,:) = 0;
%         Jy(end,:) = 0;
%     end
%     if (Jy_padding_y == 1)
%         % TRIAL
%         Jy(:,2) = Jy(:,2) + Jy(:,1);
%         Jy(:,end-1) = Jy(:,end-1) + Jy(:,end);
%         % END TRIAL
%         Jy(:,1) = 0;
%         Jy(:,end) = 0;
%     end
    % Enforce boundary conditions
%     Jy(1+Jy_padding_x,:) = 0;
%     Jy(size(Jy,1)-Jy_padding_x,:) = 0;
% catch exception
%     throw(exception);
% end