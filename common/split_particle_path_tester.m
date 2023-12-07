a_x = 0;
b_x = 1;
a_y = 0;
b_y = 1;

prev_particle_location = [.4,.2,0];
particle_location = [.6,.3,0];

[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);
assert(isequal(prev, prev_particle_location));
assert(isequal(curr, particle_location));

prev_particle_location = [.4,.2,0];
particle_location = [-.1,.3,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

% [m1,b1] = line_from_points(x1,x2,y1,y2);
% l1 = @(x_val) m1*x_val + b1;
% 
% [m2,b2] = line_from_points(x3,x4,y3,y4);
% l2 = @(x_val) m2*x_val + b2;
% 
% dx = .01;
% if (x2 < x1)
%     dx = dx*-1;
% end
% 
% xs1 = x1:dx:x2;
% xs2 = x3:dx:x4;

% plot(xs1,l1(xs1));
% hold on;
% plot(xs2,l2(xs2));
% hold off;
% drawnow;

assert(x2 == a_x);
assert(x3 == b_x);
assert(y2 == y3);

prev_particle_location = [.4,.2,0];
particle_location = [1.1,.3,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

l1 = line_from_points(x1,x2,y1,y2);
l2 = line_from_points(x3,x4,y3,y4);

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(x2 == b_x);
assert(x3 == a_x);
assert(y2 == y3);

% Moving up in a straight line (m = NaN)
prev_particle_location = [.4,.9,0];
particle_location = [.4,1.1,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

l1 = line_from_points(x1,x2,y1,y2);
l2 = line_from_points(x3,x4,y3,y4);

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(y2 == b_y);
assert(y3 == a_y);
assert(x2 == x3);

% Moving right in a straight line (m = 0)
prev_particle_location = [.4,.9,0];
particle_location = [1.3,.9,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

l1 = line_from_points(x1,x2,y1,y2);
l2 = line_from_points(x3,x4,y3,y4);

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(x2 == b_x);
assert(x3 == a_x);
assert(y2 == y3);

% Crossing b_x multiple times
prev_particle_location = [.4,.2,0];
particle_location = [2.1,.3,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

x5 = prev(3,1);
x6 = curr(3,1);
y5 = prev(3,2);
y6 = curr(3,2);

l1 = line_from_points(x1,x2,y1,y2);
% l1 = @(x_val) m1*x_val + b1;

l2 = line_from_points(x3,x4,y3,y4);
% l2 = @(x_val) m2*x_val + b2;

l3 = line_from_points(x5,x6,y5,y6);
% l3 = @(x_val) m3*x_val + b3;

dx = .01;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;
xs3 = x5:dx:x6;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
plot(xs3,l3(xs3));
hold off;
drawnow;

assert(x2 == b_x);
assert(x3 == a_x);
assert(x4 == b_x);
assert(x5 == a_x);
assert(y2 == y3);
assert(y4 == y5);

prev_particle_location = [.4,.2,0];
particle_location = [.5,-.3,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

l1 = line_from_points(x1,x2,y1,y2);
% l1 = @(x_val) m1*x_val + b1;

l2 = line_from_points(x3,x4,y3,y4);
% l2 = @(x_val) m2*x_val + b2;

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(y2 == a_y);
assert(y3 == b_y);
assert(x2 == x3);

prev_particle_location = [.4,.2,0];
particle_location = [.1,1.3,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

l1 = line_from_points(x1,x2,y1,y2);
% l1 = @(x_val) m1*x_val + b1;

l2 = line_from_points(x3,x4,y3,y4);
% l2 = @(x_val) m2*x_val + b2;

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(y2 == b_y);
assert(y3 == a_y);
assert(x2 == x3);plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(y2 == b_y);
assert(y3 == a_y);
assert(x2 == x3);

prev_particle_location = [.4,.2,0];
particle_location = [.1,2.3,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

x5 = prev(3,1);
x6 = curr(3,1);
y5 = prev(3,2);
y6 = curr(3,2);

l1 = line_from_points(x1,x2,y1,y2);
% l1 = @(x_val) m1*x_val + b1;

l2 = line_from_points(x3,x4,y3,y4);
% l2 = @(x_val) m2*x_val + b2;
l3 = line_from_points(x5,x6,y5,y6);

% l3 = @(x_val) m3*x_val + b3;

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;
xs3 = x5:dx:x6;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
plot(xs3,l3(xs3));
hold off;
drawnow;

assert(y2 == b_y);
assert(y3 == a_y);
assert(y4 == b_y);
assert(y5 == a_y);
assert(x2 == x3);
assert(x4 == x5);

prev_particle_location = [.1,.2,0];
particle_location = [-.1,-.1,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

x5 = prev(3,1);
x6 = curr(3,1);
y5 = prev(3,2);
y6 = curr(3,2);

l1 = line_from_points(x1,x2,y1,y2);
l2 = line_from_points(x3,x4,y3,y4);
l3 = line_from_points(x5,x6,y5,y6);

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;
xs3 = x5:dx:x6;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
plot(xs3,l3(xs3));
hold off;
legend("x1","x2","x3")
drawnow;

assert(y2 == y3);
assert(x4 == x5);
assert(x2 == a_x);
assert(x3 == b_x);
assert(y4 == a_y);
assert(y5 == b_y);

prev_particle_location = [.9687,.0329,0];
particle_location = [1.0031,-.0274,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

l1 = line_from_points(x1,x2,y1,y2);
% l1 = @(x_val) m1*x_val + b1;

l2 = line_from_points(x3,x4,y3,y4);
% l2 = @(x_val) m2*x_val + b2;

dx = .001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
hold off;
drawnow;

assert(y2 == a_y);
assert(y3 == b_y);
assert(x2 == x3);

prev_particle_location = [.0012,.9815,0];
particle_location = [-.0273,1.0031,0];
[prev, curr] = split_particle_path(prev_particle_location, particle_location, a_x, b_x, a_y, b_y);

x1 = prev(1,1);
x2 = curr(1,1);
y1 = prev(1,2);
y2 = curr(1,2);

x3 = prev(2,1);
x4 = curr(2,1);
y3 = prev(2,2);
y4 = curr(2,2);

x5 = prev(3,1);
x6 = curr(3,1);
y5 = prev(3,2);
y6 = curr(3,2);

l1 = line_from_points(x1,x2,y1,y2);
l2 = line_from_points(x3,x4,y3,y4);
l3 = line_from_points(x5,x6,y5,y6);

dx = .0001;
if (x2 < x1)
    dx = dx*-1;
end

xs1 = x1:dx:x2;
xs2 = x3:dx:x4;
xs3 = x5:dx:x6;

plot(xs1,l1(xs1));
hold on;
plot(xs2,l2(xs2));
plot(xs3,l3(xs3));
hold off;
legend("x1","x2","x3");
drawnow;

assert(y2 == y3);
assert(x4 == x5);
assert(x2 == a_x);
assert(x3 == b_x);
