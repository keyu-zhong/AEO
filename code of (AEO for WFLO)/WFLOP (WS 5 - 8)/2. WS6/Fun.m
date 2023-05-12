function [cost_per, P_sum, Efficiency, lp_power_accum, cr_position, xy_position] = Fun(pop, rows, cols, N, cell_width)
    theta = [0, 3 * pi / 6.0, 6 * pi / 6.0, 9 * pi / 6.0]; % 风向角
    velocity = 12.0; % 风速
    f_theta_v = [0.25; 0.25; 0.25; 0.25];% 风概率
    % fitness_val = 0;
    % P_sum = 0;
    % cost_per = 0;
    cr_position = zeros(2, N);
    xy_position = zeros(2, N);
    
    ind_position = zeros(1, N);
    ind_pos = 1;
    
    for ind = 1: rows * cols
        if pop(ind) == 1
            r_i = ceil(ind / cols);
            c_i = ceil(ind - r_i * cols);
            cr_position(1, ind_pos) = c_i;
            cr_position(2, ind_pos) = r_i;
            xy_position(1, ind_pos) = c_i * cell_width + cell_width/2;
            xy_position(2, ind_pos) = r_i * cell_width + cell_width/2;
            ind_position(ind_pos) = ind;
            ind_pos = ind_pos + 1;
        end
    end
    lp_power_accum = zeros(1,N);  % a specific layout power accumulate
    for ind_t = 1:numel(theta)
        for ind_v = 1:numel(velocity)
            trans_matrix = [cos(theta(ind_t)), -sin(theta(ind_t)); sin(theta(ind_t)), cos(theta(ind_t))];
            
            trans_xy_position = trans_matrix * xy_position;
            speed_deficiency = wake_calculate(trans_xy_position, N);
            
            actual_velocity = (1 - speed_deficiency) * velocity(ind_v);
            lp_power = layout_power(actual_velocity,N);
            
            lp_power1 = lp_power * f_theta_v(ind_t, ind_v);
            lp_power_accum = lp_power_accum + lp_power1;
            
        end
    end
    [~,sorted_index] = sort(lp_power_accum);
    % po[i, :] = ind_position[sorted_index]
    P_sum = sum(lp_power_accum);
    P_without_wake = cal_P_rate_total(theta, f_theta_v, velocity, N);
    Efficiency = P_sum  / P_without_wake;
    cost_per = cost(N, P_sum);
    P_sum = -P_sum;
    Efficiency = -Efficiency;
end

function wake_deficiency = wake_calculate(trans_xy_position, N)
    hub_height = 80.0;
    rator_radius = 77.0 * 0.5;
    surface_roughness = 0.25 * 0.001;
    entrainment_const = 0.5 / log(hub_height / surface_roughness);
    [~,sorted_index] = sort(-trans_xy_position(2, :)) ;
    wake_deficiency = zeros(1,N);
    wake_deficiency(sorted_index(1)) = 0;
    for i = 2:N
        for j = 1:i-1
            xdis = abs(trans_xy_position(1, sorted_index(i)) - trans_xy_position(1, sorted_index(j)));
            ydis = abs(trans_xy_position(2, sorted_index(i)) - trans_xy_position(2, sorted_index(j)));
            d = cal_deficiency(xdis, ydis, rator_radius,entrainment_const);
            wake_deficiency(sorted_index(i)) = wake_deficiency(sorted_index(i)) + d^2;
        end
        wake_deficiency(sorted_index(i)) = sqrt(wake_deficiency(sorted_index(i)));
    end
end


function d = cal_deficiency(dx, dy, r, ec)
    if dy == 0
        d = 0;
    else
        R = r + ec * dy;
        inter_area = cal_interaction_area(dx, dy, r, R);
        d = 2.0 / 3.0 * (r^2) / (R^2) * inter_area / (pi * r^2);
    end
end

function res = cal_interaction_area(dx, dy, r, R)
    if dx >= r + R
        res = 0;
    elseif dx >= sqrt(R^2 - r^2)
        alpha = acos((R^2 + dx^2 - r^2) / (2 * R * dx));
        beta = acos((r^2 + dx^2 - R^2) / (2 * r * dx));
        A1 = alpha * R^2;
        A2 = beta * r^2;
        A3 = R * dx * sin(alpha);
        res =  A1 + A2 - A3;
    elseif dx >= R - r
        alpha = acos((R^2 + dx^2 - r^2) / (2 * R * dx));
        beta = pi - acos((r^2 + dx^2 - R^2) / (2 * r * dx));
        A1 = alpha * R^2;
        A2 = beta * r^2;
        A3 = R * dx * sin(alpha);
        res =  pi * r^2 - (A2 + A3 - A1);
    else
        res =  pi * r^2;
    end
end

function power = layout_power(velocity, N)
    power = zeros(1,N);
    for i = 1:N
        power(i) = P_i_X(velocity(i));
    end
end

% calculate total rate power
function res = cal_P_rate_total(theta, f_theta_v, velocity, N)
    f_p = 0;
    for ind_t = 1:numel(theta)
        for ind_v = 1:numel(velocity)
            f_p = f_p + f_theta_v(ind_t, ind_v) * P_i_X(velocity(ind_v));
        end
    end
    res = N * f_p;
end

function res = P_i_X(v)
    if v < 2.0
        res = 0;
    elseif v >= 2.0 && v < 12.8
        res = 0.3 * v^3;
    elseif v >= 12.8 && v <= 18
        res = 629.1;
    else
        res = 0;
    end
end

% the total cost
function cost_per = cost(N, power_WF)
    cost_base = 1.0 * N * (2.0 / 3.0 + 1.0 / 3.0 * exp(-0.00174 * N^2));
    cost_per = 1.0 * cost_base / power_WF;
end


