index_offset = 1;    
Jx_including_ghosts = zeros(size(Jx,1)+1,size(Jx,2)+2);

x_staggered = x;
y_staggered = y_padded;

if (Jx_stagger_x ~= 0)
    x_staggered = x_padded(1:end-1)+Jx_stagger_x;
end

if (Jx_stagger_y ~= 0)
    y_staggered = y_padded(1:end-1)+Jx_stagger_y;
end
locs = particles(:,1:3);
vels = particles(:,4:6);
tags = particles(:,7);
ws = particles(:,9);
qs = particles(:,10);

lc_x = Jx_padding_x+index_offset+(locs(:,1)-Jx_stagger_x+offset_x)./dx;
lc_y = Jx_padding_y+index_offset+(locs(:,2)-Jx_stagger_y+offset_y)./dy;
is = floor(lc_x);
js = floor(lc_y);
wxs = qs.*ws.*vels(:,1);

x_nodes = x_staggered(is)';
y_nodes = y_staggered(js)';

fxs = (locs(:,1) - x_nodes)./dx;
fys = round((locs(:,2) - y_nodes)./dy);
Jx5 = accumarray([is,js],(1-fxs).*(1-fys).*wxs,size(Jx_including_ghosts));
Jx6 = accumarray([is,js+1],(1-fxs).*fys.*wxs,size(Jx_including_ghosts));
Jx7 = accumarray([is+1,js],fxs.*(1-fys).*wxs,size(Jx_including_ghosts));
Jx8 = accumarray([is+1,js+1],fxs.*fys.*wxs,size(Jx_including_ghosts));

Jx_including_ghosts = Jx5+Jx6+Jx7+Jx8;
Jx = Jx_including_ghosts(2:end,2:end-1);
Jx(1,:) = Jx(1,:) + Jx_including_ghosts(end,2:end-1);
Jx(end,:) = Jx(end,:) + Jx_including_ghosts(2,2:end-1);
% Jx(end-1,:) = Jx(end-1,:) + Jx_including_ghosts(1,2:end-1); 
% Not sure what to do with this

Jx(2:end-1,1) = Jx(2:end-1,1) + Jx_including_ghosts(3:end-1,end-1);
Jx(2:end-1,end) = Jx(2:end-1,end) + Jx_including_ghosts(3:end-1,2);

Jx(2:end-1,2) = Jx(2:end-1,2) + Jx_including_ghosts(3:end-1,end);
Jx(2:end-1,end-1) = Jx(2:end-1,end-1) + Jx_including_ghosts(3:end-1,1);
% disp(Jx);