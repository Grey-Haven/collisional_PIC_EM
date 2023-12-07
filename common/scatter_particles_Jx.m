% try
    index_offset = 1;    
    Jx_including_borders = zeros(size(Jx,1),size(Jx,2)+2);
    Jx_including_borders_vec = zeros(size(Jx,1),size(Jx,2)+2);

    x_staggered = x;
    y_staggered = y;

    if (Jx_stagger_x ~= 0)
        x_staggered = x_padded(1:end-1)+Jx_stagger_x;
    end

    if (Jx_stagger_y ~= 0)
        y_staggered = y_padded(1:end-1)+Jx_stagger_y;
    end
    
%     locs = zeros(length(particles),3);
%     vels = zeros(length(particles),3);
%     qs = zeros(length(particles),1);
%     ws = zeros(length(particles),1);
    
%     for p = 1:length(particles)
%         particle = particles{p};
%         locs(p,:) = particle.location;
%         vels(p,:) = particle.velocity;
%         qs(p) = particle.charge;
%         ws(p) = weightMap(particle.tag);
%     end

    locs = particles(:,1:3);
    vels = particles(:,4:6);
%     tags = particles(:,7);
%     ms = particles(:,8);
    ws = particles(:,9);
    qs = particles(:,10);
%     ws = weightMap(char(tags));

    lc_x = Jx_padding_x+index_offset+(locs(:,1)-Jx_stagger_x+offset_x)./dx;
    lc_y = Jx_padding_y+index_offset+(locs(:,2)-Jx_stagger_y+offset_y)./dy;
    is = floor(lc_x);
    js = floor(lc_y);
    wxs = qs.*ws.*vels(:,1);
    
    x_nodes = x_staggered(is)';
    y_nodes = y_staggered(js)';

    fxs = (locs(:,1) - x_nodes)./dx;
    fys = (locs(:,2) - y_nodes)./dy;
    
%     idxs = sub2ind(size(Jx_including_borders_vec),is,js);
%     Jx_including_borders_vec(idxs) = Jx_including_borders_vec(idxs) + (1-fxs).*(1-fys).*wxs;
%     for p = 1:length(particles)
%         i = is(p);
%         j = js(p);
%         w = wxs(p);
%         fx = fxs(p);
%         fy = fys(p);
%         Jx_including_borders_vec(i,j) = Jx_including_borders_vec(i,j) + (1-fx).*(1-fy).*w;
%     end
    
%     Jx1 = full(sparse(is,js,(1-fxs).*(1-fys).*wxs,size(Jx,1),size(Jx,2)));
%     Jx2 = full(sparse(is,js+1,(1-fxs).*fys.*wxs,size(Jx,1),size(Jx,2)));
%     Jx3 = full(sparse(is+1,js,fxs.*(1-fys).*wxs,size(Jx,1),size(Jx,2)));
%     Jx4 = full(sparse(is+1,js+1,fxs.*fys.*wxs,size(Jx,1),size(Jx,2)));
    
    Jx5 = accumarray([is,js],(1-fxs).*(1-fys).*wxs,size(Jx_including_borders));
    Jx6 = accumarray([is,js+1],(1-fxs).*fys.*wxs,size(Jx_including_borders));
    Jx7 = accumarray([is+1,js],fxs.*(1-fys).*wxs,size(Jx_including_borders));
    Jx8 = accumarray([is+1,js+1],fxs.*fys.*wxs,size(Jx_including_borders));
    
%     Jx_including_borders_vec1 = Jx1+Jx2+Jx3+Jx4;
    Jx_including_borders_vec2 = Jx5+Jx6+Jx7+Jx8;
    Jx_including_borders = Jx_including_borders_vec2;
        
%     Jx_including_borders_vec(is,js) = Jx_including_borders_vec(is,js) + (1-fxs).*(1-fys).*wxs;
%     Jx_including_borders_vec(is,js+1) = Jx_including_borders_vec(is,js+1) + (1-fxs).*fys.*wxs;
%     Jx_including_borders_vec(is+1,js) = Jx_including_borders_vec(is+1,js) + fxs.*(1-fys).*wxs;
%     Jx_including_borders_vec(is+1,js+1) = Jx_including_borders_vec(is+1,js+1) + fxs.*fys.*wxs;
    
%     for p = 1:size(particles,1)
%        particle = particles(p,:);
%         
%         locs = particle(1:3);
%         vels = particle(4:6);
% %         tags = particle(7);
%         weight = particle(9);
%         q = particle(10);
% 
%         lc_x = Jx_padding_x+index_offset+(locs(1)-Jx_stagger_x+offset_x)/dx;
%         lc_y = Jx_padding_y+index_offset+(locs(2)-Jx_stagger_y+offset_y)/dy;
%                 
%         i = floor(lc_x);
%         j = floor(lc_y);
% 
% %         q = particle.charge;
%         vx = vels(1);
% %         weight = weightMap(particle.tag);
%         wx = q*vx*weight;
% 
%         x_staggered = x;
%         y_staggered = y;
% 
%         if (Jx_stagger_x ~= 0)
%             x_staggered = x_padded(1:end-1)+Jx_stagger_x;
%         end
% 
%         if (Jx_stagger_y ~= 0)
%             y_staggered = y_padded(1:end-1)+Jx_stagger_y;
%         end
%         
%         x_node = x_staggered(i); % - Jx_stagger_x;
%         y_node = y_staggered(j); % - Jx_stagger_y;
% 
%         fx = (locs(1) - x_node)/dx;
%         fy = (locs(2) - y_node)/dy;
%         
%         if (fx < 0 || fy < 0 || fx > 1 || fy > 1)
%             throw(exception);
%         end
%         
%         Jx_including_borders(i,j) = Jx_including_borders(i,j) + (1-fx)*(1-fy)*wx;
%         Jx_including_borders(i,j+1) = Jx_including_borders(i,j+1) + (1-fx)*fy*wx;
%         Jx_including_borders(i+1,j) = Jx_including_borders(i+1,j) + fx*(1-fy)*wx;
%         Jx_including_borders(i+1,j+1) = Jx_including_borders(i+1,j+1) + fx*fy*wx;
%     end
%     eps = 1e-15;
%     if ~(abs(norm(Jx_including_borders-Jx_including_borders_vec2)) < eps)
%         disp('foo');
%     end
%     assert(isequal(Jx_including_borders,Jx_including_borders_vec2));
    Jx_including_borders(2,:) = Jx_including_borders(2,:) + Jx_including_borders(1,:);
    Jx_including_borders(end-1,:) = Jx_including_borders(end-1,:) + Jx_including_borders(end,:);
    Jx(:,:) = Jx_including_borders(:,2:end-1); % cut out boundaries (perfectly conducting)
%     Jx(2,:) = Jx(2,:) + Jx(1,:);
%     Jx(end-1,:) = Jx(end-1,:) + Jx(end,:);
%     if (Jx_padding_x == 1)
% %         TRIAL
%         Jx(2,:) = Jx(2,:) + Jx(1,:);
%         Jx(end-1,:) = Jx(end-1,:) + Jx(end,:);
%         % END TRIAL
%         Jx(1,:) = 0;
%         Jx(end,:) = 0;
%     end
%     if (Jx_padding_y == 1)
%         % TRIAL
%         Jx(:,2) = Jx(:,2) + Jx(:,1);
%         Jx(:,end-1) = Jx(:,end-1) + Jx(:,end);
%         % END TRIAL
%         Jx(:,1) = 0;
%         Jx(:,size(Jx,2)) = 0;
%     end
    % Enforce boundary conditions
%     Jx(:,1+Jx_padding_y) = 0;
%     Jx(:,size(Jx,2)-Jx_padding_y) = 0;
% catch exception
%     throw(exception);
% end