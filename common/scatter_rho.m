function rho = scatter_rho(locs,qs,ws,x,y)

    index_offset = 1;
%     pad_x = 1;
%     pad_y = 1;

    % x_padded = x(1)-dx:dx:x(end)+dx;
    % y_padded = y(1)-dy:dy:y(end)+dy;
    
    offset_x = -x(1);
    offset_y = -y(1);
    
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
    rho = zeros(length(x),length(y));

    % if (Hz_stagger_x ~= 0)
    %     x_staggered = x_padded(1:end-1)+Hz_stagger_x;
    % end
    % 
    % if (Hz_stagger_y ~= 0)
    %     y_staggered = y_padded(1:end-1)+Hz_stagger_y;
    % end

%     locs = particles(:,1:3);
%     vels = particles(:,4:6);
%     tags = particles(:,7);
%     ms = particles(:,8);
%     ws = particles(:,9);
%     qs = particles(:,10);
    %     ws = weightMap(char(tags));

    lc_x = index_offset+(locs(:,1)+offset_x)./dx;
    lc_y = index_offset+(locs(:,2)+offset_y)./dy;
    is = floor(lc_x);
    js = floor(lc_y);
    wzs = qs.*ws;

    x_nodes = x(is)';
    y_nodes = y(js)';

    fxs = (locs(:,1) - x_nodes)/dx;
    fys = (locs(:,2) - y_nodes)/dy;

    assert(all(0 <= fxs & fxs <= 1))

    rho5 = accumarray([is,js],(1-fxs).*(1-fys).*wzs,[size(rho,1)+1,size(rho,2)]);
    rho6 = accumarray([is,js+1],(1-fxs).*fys.*wzs,[size(rho,1)+1,size(rho,2)]);
    rho7 = accumarray([is+1,js],fxs.*(1-fys).*wzs,[size(rho,1)+1,size(rho,2)]);
    rho8 = accumarray([is+1,js+1],fxs.*fys.*wzs,[size(rho,1)+1,size(rho,2)]);

    rho = (rho5(1:end-1,:)+rho6(1:end-1,:)+rho7(1:end-1,:)+rho8(1:end-1,:)) ./ (dx*dy);
%     rho = rho5(1:end-1,:)+rho6(1:end-1,:)+rho7(1:end-1,:)+rho8(1:end-1,:);

end