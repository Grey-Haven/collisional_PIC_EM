index_offset = 1;
Jy_including_ghosts = zeros(size(Jy,1)+2,size(Jy,2)+1);

x_staggered = x_padded;
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

lc_x = Jy_padding_x+index_offset+(locs(:,1)-Jy_stagger_x+offset_x)./dx;
lc_y = Jy_padding_y+index_offset+(locs(:,2)-Jy_stagger_y+offset_y)./dy;
is = floor(lc_x);
js = floor(lc_y);
wys = qs.*ws.*vels(:,2);

x_nodes = x_staggered(is)';
y_nodes = y_staggered(js)';

fxs = (locs(:,1) - x_nodes)./dx;
fys = (locs(:,2) - y_nodes)./dy;

Jy5 = accumarray([is,js],(1-fxs).*(1-fys).*wys,size(Jy_including_ghosts));
Jy6 = accumarray([is,js+1],(1-fxs).*fys.*wys,size(Jy_including_ghosts));
Jy7 = accumarray([is+1,js],fxs.*(1-fys).*wys,size(Jy_including_ghosts));
Jy8 = accumarray([is+1,js+1],fxs.*fys.*wys,size(Jy_including_ghosts));

Jy_including_ghosts = Jy5+Jy6+Jy7+Jy8;
Jy = Jy_including_ghosts(2:end-1,2:end);

Jy(:,1) = Jy(:,1) + Jy_including_ghosts(2:end-1,end);
Jy(:,end) = Jy(:,end) + Jy_including_ghosts(2:end-1,2);

Jy(1,2:end-1) = Jy(1,2:end-1) + Jy_including_ghosts(end-1,3:end-1);
Jy(end,2:end-1) = Jy(end,2:end-1) + Jy_including_ghosts(2,3:end-1);

Jy(2,2:end-1) = Jy(2,2:end-1) + Jy_including_ghosts(end,3:end-1);
Jy(end-1,2:end-1) = Jy(end-1,2:end-1) + Jy_including_ghosts(1,3:end-1);
% disp(Jy);