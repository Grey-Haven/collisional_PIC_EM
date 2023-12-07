function [x_prev_trim, y_prev_trim, x_curr_trim, y_curr_trim, dt_trim] = trim_particle_path(x_prev,y_prev,x_curr,y_curr,x,y,dt)

%     x_prev_trim = x_prev;
%     y_prev_trim = y_prev;
%     x_curr_trim = x_curr;
%     y_curr_trim = y_curr;
    
%     dt_trim = dt;

    [line,line_inv] = line_from_points(x_prev,x_curr,y_prev,y_curr);

    % Assumption: the particle will only cross two of four barriers, and
    % these must be connected (ie exit via bottom left, enter via top
    % right, etc). ie it will never enter left and exit right.
    
    % Possible cases:    
    % particle exits/enters right
    % particle exits/enters left
    % particle exits/enters up
    % particle exits/enters down
    % particle exits/enters first right then up
    % particle exits/enters first right then down
    % particle exits/enters first up then right
    % particle exits/enters first down then right
    % particle exits/enters first left then up
    % particle exits/enters first left then down
    % particle exits/enters first up then left
    % particle exits/enters first down then left
    % (The last eight cases should absorb the (literal) corner cases)
    
    x_prev_trim = max(x(1), min(x(end), x_prev));
    x_curr_trim = max(x(1), min(x(end), x_curr));
    y_prev_trim0 = line(x_prev_trim);
    y_curr_trim0 = line(x_curr_trim);
    % Case when particle is going straight up/down
    if (~isnan(y_prev_trim0))
        y_prev_trim = y_prev_trim0;
    else
        y_prev_trim = y_prev;
    end
    if (~isnan(y_curr_trim0))
        y_curr_trim = y_curr_trim0;
    else
        y_curr_trim = y_curr;
    end
    
    y_prev_trim0 = max(y(1), min(y(end), y_prev_trim));
    y_curr_trim0 = max(y(1), min(y(end), y_curr_trim));
    x_prev_trim0 = line_inv(y_prev_trim0);
    x_curr_trim0 = line_inv(y_curr_trim0);
    
    if (y_prev_trim < y(1) || y_prev_trim > y(end))
        y_prev_trim = y_prev_trim0;
        x_prev_trim = x_prev_trim0;
    end
    if (y_curr_trim < y(1) || y_curr_trim > y(end))
        y_curr_trim = y_curr_trim0;
        x_curr_trim = x_curr_trim0;
    end
    
    % Case when particle is going straight left/right
%     if (~isnan(x_prev_trim0))
%         x_prev_trim = x_prev_trim0;
%     end
%     if (~isnan(x_curr_trim0))
%         x_curr_trim = x_curr_trim0;
%     end
    
    del_x = abs(x_curr - x_prev);
    del_y = abs(y_curr - y_prev);
    del_s = sqrt(del_x^2 + del_y^2);
    
    del_x_trim = abs(x_curr_trim - x_prev_trim);
    del_y_trim = abs(y_curr_trim - y_prev_trim);
    del_s_trim = sqrt(del_x_trim^2 + del_y_trim^2);
    
    frac = del_s_trim/del_s;
    dt_trim = frac*dt;
end