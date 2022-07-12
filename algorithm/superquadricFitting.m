function [x, cost, residual] = superquadricFitting(point, para, x0)

iter_max = para.iterMax;
iter_min = para.iterMin;
tolerance = para.tolerance;
relative_tolerance = para.relative_tolerance;
adaptive_upper = para.adaptive_upper;
taper = 1;
% optimization settings
options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'Display', 'off', 'MaxIterations', ...
    2, 'SpecifyObjectiveGradient', true);

%% Initialization
% set lower and upper bounds for the superquadrics
upper = 4 * max(max(abs(point)));
lb = [0.0 0.0 0.00001 0.00001 0.00001 -2*pi -2*pi -2*pi -ones(1, 3) * upper + x0(9:11) -taper -taper];
ub = [2.0 2.0 ones(1, 3) * upper  2*pi 2*pi 2*pi ones(1, 3) * upper + x0(9:11) taper taper];

%% optimization
% initialize parameters
x = x0;
cost = inf;
switched = 0;

% optimization
for iter = 1 : iter_max

    if adaptive_upper == 1
        R_current = eul2rotm(x(6 : 8));
        point_rot_current = R_current' * point - R_current' * x(9 : 11)';
%         ub_a = s * [max(point_rot_current(1, pIdx))-min(point_rot_current(1, pIdx)), max(point_rot_current(2, pIdx))-min(point_rot_current(2, pIdx)),max(point_rot_current(3, pIdx))-min(point_rot_current(3,pIdx))];
        ub_a = [0.5 * (max(point_rot_current(1, :)) - min(point_rot_current(1, :))), ...
                0.5 * (max(point_rot_current(2, :)) - min(point_rot_current(2, :))), ...
                0.5 * (max(point_rot_current(3, :)) - min(point_rot_current(3, :)))];
        ub_a(ub_a < 0.00001) = 0.00001;
        
        ub = [2.0 2.0 ub_a  2*pi 2*pi 2*pi abs(ub_a*R_current') + x(9 : 11) taper taper];
        lb = [0.0 0.0 0.00001 0.00001 0.00001 -2*pi -2*pi -2*pi -abs(ub_a*R_current') + x(9 : 11) -taper -taper];
    end
    
    cost_func = @(x) weighted_dist(x, point, ones(1,size(point,2)));
    [x_n, cost_n, residual_n] = lsqnonlin(cost_func, x, lb, ub, options);
    % evaluate relative cost decrease
    relative_cost = (cost - cost_n) / cost_n;

    if (cost_n < tolerance && iter > 1) || (relative_cost < relative_tolerance && switched >= para.max_switch && iter > iter_min) % >1 >5
        cost = cost_n;
        x = x_n;
        residual = residual_n;
        break
    end

    if relative_cost < relative_tolerance && iter ~= 1 % set different tolerance for switch and termination
        % activate switching algorithm to avoid local minimum
        switch_success = 0;
        % case1 - axis-mismatch similarity
        axis_0 = eul2rotm(x(6 : 8));
        axis_1 = circshift(axis_0, [0, 2]);
        axis_2 = circshift(axis_0, [0, 1]);
        eul_1 = rotm2eul(axis_1);
        eul_2 = rotm2eul(axis_2);
        x_axis = [x(2), x(1), x(4), x(5), x(3), eul_1, x(9 : 13); ...
                  x(2), x(1), x(5), x(3), x(4), eul_2, x(9 : 13)];
        % case2 - duality similarity
        scale_ratio = circshift(x(3 : 5), 2) ./ x(3 : 5);
        scale_idx = find(and(scale_ratio > 0.9, scale_ratio < 1.1)); %0.5 1.5
        x_rot = zeros(size(scale_idx, 2), 13);
        rot_idx = 1;
        if ismember(1, scale_idx)
            eul_rot = rotm2eul(axis_0 * rotz(45));
            if x(2) <= 1
                x_rot(rot_idx, :) = [x(1), 2 - x(2), ((1 - sqrt(2)) * x(2) + sqrt(2)) * min(x(3), x(4)) * ones(1, 2), x(5), eul_rot, x(9 : 13)];
            else
                x_rot(rot_idx, :) = [x(1), 2 - x(2), ((sqrt(2)/2 - 1) * x(2) + 2 - sqrt(2)/2) * min(x(3), x(4)) * ones(1, 2), x(5), eul_rot, x(9 : 13)];
            end            
            rot_idx = rot_idx + 1;
        end
        if ismember(2, scale_idx)
            eul_rot = rotm2eul(axis_1 * rotz(45));
            if x(1) <= 1
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(4), x(5)) * ones(1, 2), x(3), eul_rot, x(9 : 13)];
            else
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((sqrt(2)/2 - 1) * x(1) + 2 - sqrt(2)/2) * min(x(4), x(5)) * ones(1, 2), x(3), eul_rot, x(9 : 13)];
            end    
            rot_idx = rot_idx + 1;
        end
        if ismember(3, scale_idx)
            eul_rot = rotm2eul(axis_2 * rotz(45));       
            if x(1) <= 1
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(5), x(3)) * ones(1, 2), x(4), eul_rot, x(9 : 13)];
            else
                x_rot(rot_idx, :) = [x(2), 2 - x(1), ((sqrt(2)/2 - 1) * x(1) + 2 - sqrt(2)/2) * min(x(5), x(3)) * ones(1, 2), x(4), eul_rot, x(9 : 13)];
            end    
        end

        % generate candidate configuration list with cost
        x_candidate = [x_axis; x_rot];
        cost_candidate = weighted_cost_switch(x_candidate, point, ones(1,size(point,2)));

        idx_nan = find(and(~isnan(cost_candidate), ~isinf(cost_candidate)));
        cost_candidate = cost_candidate(idx_nan);
        x_candidate = x_candidate(idx_nan, :);
        [~, idx] = sort(cost_candidate);

        for i_candidate = 1 : size(idx, 1)
                % adaptive upper bound
            if adaptive_upper == 1
                R_current = eul2rotm(x_candidate(idx(i_candidate), 6 : 8));
                point_rot_current = R_current' * point - R_current' * x_candidate(idx(i_candidate), 9 : 11)';
%                 ub_a = s * [max(point_rot_current(1, pIdx))-min(point_rot_current(1, pIdx)), max(point_rot_current(2, pIdx))-min(point_rot_current(2, pIdx)),max(point_rot_current(3, pIdx))-min(point_rot_current(3,pIdx))];
                ub_a = [0.5 * (max(point_rot_current(1, :)) - min(point_rot_current(1, :))), ...
                        0.5 * (max(point_rot_current(2, :)) - min(point_rot_current(2, :))), ...
                        0.5 * (max(point_rot_current(3, :)) - min(point_rot_current(3, :)))];
                ub_a(ub_a < 0.00001) = 0.00001;
                ub = [2.0 2.0 ub_a  2*pi 2*pi 2*pi abs(ub_a*R_current') + x(9:11) taper taper];
                lb = [0.0 0.0 0.00001 0.00001 0.00001 -2*pi -2*pi -2*pi -abs(ub_a*R_current') + x(9:11) -taper -taper];
            end
            [x_switch, cost_switch, residual_switch] = lsqnonlin(cost_func, x_candidate(idx(i_candidate), :), lb, ub, options);
            if cost_switch < min(cost_n, cost)
                x = x_switch;
                cost = cost_switch;
                residual = residual_switch;
                % update sigma
                switch_success = 1;
                break
            end
        end
        if switch_success == 0
            cost = cost_n;
            x = x_n;
            residual = residual_n;
        end
        switched = switched + 1;
    else
        cost = cost_n;
        x = x_n;
        residual = residual_n;
    end

end

end


%% Functions

% ------------------weighed distance function -----------------------------
function [value,J] = weighted_dist(para, X, p)
    [D,J] = radial_distance_with_gradient(X, para);
%         [D] = radial_distance(X, para);
    value = p .^ (1 / 2) .* D;
    value = value';

end

% ------------------weighed distance function (switch)---------------------
function [value] = weighted_cost_switch(para, X, p)
    value = zeros(size(para, 1), 1);
    for i = 1 : size(para, 1)
        value(i, 1) = sum(p .* (distance(X, para(i, :)) .^ 2), 2);
    end
end
