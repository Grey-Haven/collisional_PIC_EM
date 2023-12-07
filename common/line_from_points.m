function [line,line_inv] = line_from_points(x1,x2,y1,y2)
    m = (y2 - y1)./(x2 - x1);
    b = y2 - m.*x2;

    if (abs(m) == Inf)
        line = @(x_val) NaN;
        line_inv = @(y_val) x1;
%     elseif (isnan(m))
%         line = @(x_val) y1;
%         line_inv = @(y_val) x1;
    else
        line = @(x_val) m.*x_val+b;
        line_inv = @(y_val) (y_val-b)./m;
    end
end

