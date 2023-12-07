function [prev_particle_locations,particle_locations] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y)
    x1 = prev_particle_location(1);
    x2 = particle_location(1);
    y1 = prev_particle_location(2);
    y2 = particle_location(2);
    z1 = prev_particle_location(3);
    z2 = particle_location(3);

    [line, line_inv] = line_from_points(x1,x2,y1,y2);
    % We have (x_s, a_y) and (a_x, y_s)

    % One tricky thing, if a particle crosses the x boundary, it's its y
    % value that has exceeded the bounds of (a_y, b_y). Vice versa for y
    % boundary and x value.

    % assumption, the particle will always begin in the box
    if (x2 < a_x && y2 < a_y)
        x_s = line_inv(a_y);
        % Find out which is crossed first, a_x or a_y
        if (x_s > a_x) % particle crosses x-axis (a_y) first
            [prev_particle_locations,particle_locations] = crossed_south(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        else % particle crosses y-axis (a_x) first
            [prev_particle_locations,particle_locations] = crossed_west(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        end
    elseif (x2 > b_x && y2 < a_y)
        x_s = line_inv(a_y);
        % Find out which is crossed first, a_x or b_y
        if (x_s < b_x) % particle crosses x-axis (a_y) first
            [prev_particle_locations,particle_locations] = crossed_south(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        else % particle crosses y-axis (b_x) first
            [prev_particle_locations,particle_locations] = crossed_east(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        end
    elseif (x2 < a_x && y2 > b_y)
        x_s = line_inv(b_y);
        % Find out which is crossed first, b_x or a_y
        if (x_s > a_x) % particle crosses x-axis (b_y) first
            [prev_particle_locations,particle_locations] = crossed_north(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        else % particle crosses y-axis (a_x) first
            [prev_particle_locations,particle_locations] = crossed_west(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        end
    elseif (x2 > b_x && y2 > b_y)
        x_s = line_inv(b_y);
        % Find out which is crossed first, b_x or a_y
        if (x_s < b_x) % particle crosses x-axis (b_y) first
            [prev_particle_locations,particle_locations] = crossed_north(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        else % particle crosses y-axis (b_x) first
            [prev_particle_locations,particle_locations] = crossed_east(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
        end
    elseif x2 < a_x
%         newStartLoc = [b_x, line(a_x), 0];
%         newEndLoc = [x2 + (b_x - a_x), y2, 0];
%         [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
%         prev_particle_locations = [[x1,y1,z1];...
%                                     prev];
%         particle_locations = [[a_x,line(a_x),z2];...
%                                curr];
        [prev_particle_locations, particle_locations] = crossed_west(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
    elseif x2 > b_x
%         newStartLoc = [a_x, line(b_x), 0];
%         newEndLoc = [a_x + (x2 - b_x), y2, 0];
%         [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
%         prev_particle_locations = [[x1,y1,z1];prev];
%         particle_locations = [[b_x,line(b_x),z2];curr];
        
        [prev_particle_locations, particle_locations] = crossed_east(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
    elseif y2 < a_y
%         newStartLoc = [line_inv(a_y), b_y, 0];
%         newEndLoc = [x2, y2 + (b_y - a_y), 0];
%         [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
%         prev_particle_locations = [[x1,y1,z1];prev];
%         particle_locations = [[line_inv(a_y),a_y,z2];curr];
        [prev_particle_locations, particle_locations] = crossed_south(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
    elseif y2 > b_y
%         newStartLoc = [line_inv(b_y), a_y, 0];
%         newEndLoc = [x2, a_y + (y2 - b_y), 0];
%         [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
%         prev_particle_locations = [[x1,y1,z1];prev];
%         particle_locations = [[line_inv(b_y),b_y,z2];curr];
        [prev_particle_locations, particle_locations] = crossed_north(a_x,b_x,a_y,b_y,prev_particle_location,particle_location);
    else
        prev_particle_locations = [[x1,y1,z2]];
        particle_locations = [[x2,y2,z2]];        
    end
end

function [prev_particle_locations,particle_locations] = crossed_south(a_x,b_x,a_y,b_y,prev_particle_location,particle_location)
    [x1,x2,y1,y2,z1,z2] = unpack_particle_points(prev_particle_location,particle_location);
    [line, line_inv] = line_from_points(x1,x2,y1,y2);

    newStartLoc = [line_inv(a_y), b_y-eps, 0];
    newEndLoc = [x2, y2 + (b_y - a_y), 0];
    [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
    prev_particle_locations = [[x1,y1,z1];prev];
    particle_locations = [[line_inv(a_y),a_y,z2];curr];
end

function [prev_particle_locations,particle_locations] = crossed_west(a_x,b_x,a_y,b_y,prev_particle_location,particle_location)
        [x1,x2,y1,y2,z1,z2] = unpack_particle_points(prev_particle_location,particle_location);
        [line, line_inv] = line_from_points(x1,x2,y1,y2);
        
        newStartLoc = [b_x-eps, line(a_x), 0];
        newEndLoc = [x2 + (b_x - a_x), y2, 0];
        [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
        prev_particle_locations = [[x1,y1,z1];prev];
        particle_locations = [[a_x,line(a_x),z2];curr];
end

function [prev_particle_locations,particle_locations] = crossed_east(a_x,b_x,a_y,b_y,prev_particle_location,particle_location)
    [x1,x2,y1,y2,z1,z2] = unpack_particle_points(prev_particle_location,particle_location);
    [line, line_inv] = line_from_points(x1,x2,y1,y2);
    
    newStartLoc = [a_x, line(b_x), 0];
    newEndLoc = [a_x + (x2 - b_x), y2, 0];
    [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
    prev_particle_locations = [[x1,y1,z1];prev];
    particle_locations = [[b_x-eps,line(b_x),z2];curr];
end

function [prev_particle_locations,particle_locations] = crossed_north(a_x,b_x,a_y,b_y,prev_particle_location,particle_location)
    [x1,x2,y1,y2,z1,z2] = unpack_particle_points(prev_particle_location,particle_location);
    [line, line_inv] = line_from_points(x1,x2,y1,y2);
    
    newStartLoc = [line_inv(b_y), a_y, 0];
    newEndLoc = [x2, a_y + (y2 - b_y), 0];
    [prev, curr] = split_particle_path(newStartLoc, newEndLoc, a_x, b_x, a_y, b_y);
    prev_particle_locations = [[x1,y1,z1];prev];
    particle_locations = [[line_inv(b_y),b_y-eps,z2];curr];
end

function [x1,x2,y1,y2,z1,z2] = unpack_particle_points(prev_particle_location,particle_location)
    x1 = prev_particle_location(1);
    x2 = particle_location(1);
    y1 = prev_particle_location(2);
    y2 = particle_location(2);
    z1 = prev_particle_location(3);
    z2 = particle_location(3);
end
