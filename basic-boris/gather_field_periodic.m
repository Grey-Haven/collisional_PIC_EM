function F_p = gather_field_periodic(F,loc,x,y,z,dx,dy,dz,...
                                     stagger_x,stagger_y,stagger_z,...
                                     pad_x,pad_y,pad_z,...
                                     a_x,a_y,a_z,b_x,b_y,b_z)
                        
                        
    index_offset = 1;

    lc_x = index_offset+pad_x+(loc(:,1)-a_x+stagger_x)/dx;
    lc_y = index_offset+pad_y+(loc(:,2)-a_y+stagger_y)/dy;
    lc_z = index_offset+pad_z+(loc(:,3)-a_z+stagger_z)/dz;

    is = mod(floor(lc_x)-1,length(x))+1;
    js = mod(floor(lc_y)-1,length(y))+1;
    ks = mod(floor(lc_z)-1,length(z))+1;

    x_node = x(is)';
    y_node = y(js)';
    z_node = z(ks)';

    d1 = (loc(:,1) - x_node)./dx;
    d2 = (loc(:,2) - y_node)./dy;
    d3 = (loc(:,3) - z_node)./dz;
    
%     if (d1 < 0 || d2 < 0 || d3 < 0 || d1 > 1 || d2 > 1 || d3 > 1)
%         throw(exception);
%     end
    
%     try

    idx1 = sub2ind(size(F),is,js,ks);
    idx2 = sub2ind(size(F),is+1,js,ks);
    idx3 = sub2ind(size(F),is,js+1,ks);
    idx4 = sub2ind(size(F),is+1,js+1,ks);
    idx5 = sub2ind(size(F),is,js,ks+1);
    idx6 = sub2ind(size(F),is+1,js,ks+1);
    idx7 = sub2ind(size(F),is,js+1,ks+1);
    idx8 = sub2ind(size(F),is+1,js+1,ks+1);    

    F1 = F(idx1).*(1-d1).*(1-d2).*(1-d3);
    F2 = F(idx2).*d1.*(1-d2).*(1-d3);
    F3 = F(idx3).*(1-d1).*d2.*(1-d3);
    F4 = F(idx4).*d1.*d2.*(1-d3);
    F5 = F(idx5).*(1-d1).*d2.*d3;
    F6 = F(idx6).*d1.*(1-d2).*d3;
    F7 = F(idx7).*(1-d1).*d2.*d3;
    F8 = F(idx8).*d1.*d2.*d3;
    
    F_p = F1+F2+F3+F4+F5+F6+F7+F8;

%     F_p = F(i,j,k).*(1-d1).*(1-d2).*(1-d3);  %contribution from (i,j,k)
%     F_p = F_p + F(i+1,j,k).*d1.*(1-d2).*(1-d3);  %(i+1,j,k)
%     F_p = F_p + F(i,j+1,k).*(1-d1).*d2.*(1-d3);  %(i,j+1,k)
%     F_p = F_p + F(i+1,j+1,k).*d1.*d2.*(1-d3);  %(i+1,j+1,k)
% 
%     F_p = F_p + F(i+1,j,k+1).*(1-d1).*d2.*d3;
%     F_p = F_p + F(i+1,j,k+1).*d1.*(1-d2).*d3;  %(i+1,j,k)
%     F_p = F_p + F(i,j+1,k+1).*(1-d1).*d2.*d3;  %(i,j+1,k)
%     F_p = F_p + F(i+1,j+1,k+1).*d1.*d2.*d3;  %(i+1,j+1,k)
%     catch exception
%         throw(exception);
%     end
    for i = 1:length(F_p)
        if isnan(F_p(i))
            disp('foo');
        end
    end
end