function T_e = get_heat(vx, vy)
    avg_vx = get_mean(vx);
    avg_vy = get_mean(vy);
    var_vx = get_var(vx,avg_vx);
    var_vy = get_var(vy,avg_vy);
    avg_v = (var_vx + var_vy)/2;
    T_e = avg_v;
end

function mu = get_mean(vs)
    mu = sum(vs)/length(vs);
end

function sig = get_var(vs, mu)
    sig = sqrt(sum((vs - mu).^2)/length(vs));
end