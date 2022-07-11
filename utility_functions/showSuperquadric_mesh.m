function [] = showSuperquadric_mesh(x, varargin)

R = eul2rotm(x(6 : 8));
t = x(9 : 11);

x(1:2) = max(x(1:2),0.05);

% with tapering or not
taper = false;
color = 'r';
ViewAxis = [0 0];
CamRoll = 0;
ShowAxis = 0;
arclength = 0.1;
FaceAlpha = 1;
FaceLighting = 'flat';
lighting = false;
scale = 1;

for k = 1 : size(varargin, 2)
    if strcmp(varargin{k}, 'Taper')
        taper = varargin{k + 1};
    end
    if strcmp(varargin{k}, 'Color')
        color = varargin{k + 1};
    end
    if strcmp(varargin{k}, 'ViewAxis')
        ViewAxis = varargin{k + 1};
    end
    if strcmp(varargin{k}, 'CamRoll')
        CamRoll = varargin{k + 1};
    end
    if strcmp(varargin{k}, 'ShowAxis')
        ShowAxis= varargin{k + 1};
    end
    if strcmp(varargin{k}, 'Arclength')
        arclength= varargin{k + 1};
    end
    if strcmp(varargin{k}, 'FaceAlpha')
        FaceAlpha= varargin{k + 1};
    end
    if strcmp(varargin{k}, 'FaceLighting')
        FaceLighting= varargin{k + 1};
    end
    if strcmp(varargin{k}, 'Light')
        lighting= varargin{k + 1};
    end
    if strcmp(varargin{k}, 'Scale')
        scale = varargin{k+1};
    end
end

% validate dimensionality
if taper == true
    if size(x, 2) ~= 13
        error('Input parameters should have dimension (:, 13) for taperred SQ.')
    end
else
    if size(x, 2) ~= 11
        error('Input parameters should have dimension (:, 11) for taperred SQ.')
    end
end

[point_eta] = uniformSampledSuperellipse(x(1), [1, x(5)], arclength);
[point_omega] = uniformSampledSuperellipse(x(2), [x(3), x(4)], arclength);

x_mesh = ones(size(point_omega, 2), size(point_eta, 2));
y_mesh = ones(size(point_omega, 2), size(point_eta, 2));
z_mesh = ones(size(point_omega, 2), size(point_eta, 2));

for m = 1 : size(point_omega, 2)
    for n = 1 : size(point_eta, 2)
        point_temp = [point_omega(:, m) * point_eta(1, n); point_eta(2, n)];
        
        fx = x(12) * point_temp(3) / x(5) + 1;
        fy = x(13) * point_temp(3) / x(5) + 1;
        fz = 1;
        
        point_temp(1) = point_temp(1) * fx;
        point_temp(2) = point_temp(2) * fy;
        point_temp(3) = point_temp(3) * fz;
        
        point_temp = R * point_temp + t';
        
        x_mesh(m, n) = point_temp(1);
        y_mesh(m, n) = point_temp(2);
        z_mesh(m, n) = point_temp(3);
    end
end

x_mesh = x_mesh / scale;
y_mesh = y_mesh / scale;
z_mesh = z_mesh / scale;

mesh(x_mesh, y_mesh, z_mesh, 'FaceAlpha', FaceAlpha, 'facecolor', color, ...
    'LineStyle', 'none', 'FaceLighting', FaceLighting)
if lighting == 1
    light
    material dull
end

axis equal
view(ViewAxis)
camroll(CamRoll)

if ShowAxis == 0
    axis off
end

% hold off
% ---------------------------------utility functions ----------------------
    function [point, theta] = uniformSampledSuperellipse(epsilon, scale, arclength)
        threshold = 1e-2;
        num_limit = 10000;
        theta = zeros(1, num_limit);
        theta(1) = 0;
        
        for i = 2 : num_limit
            dt = dtheta(theta(i - 1), arclength, threshold, scale, epsilon);
            theta_temp = theta(i - 1) + dt;
            
            if theta_temp > pi/4
                break
            else
                if i < num_limit
                    theta(i) = theta_temp;
                else
                    error(['The number of the sampled points exceeds the limit of ', ...
                        num2str(num_limit * 4),...
                        '. Please increase the arclength or raise the limit'])
                end
            end
        end
        critical = i;
        
        for j = critical + 1 : num_limit
            dt = dtheta(theta(j - 1), arclength, threshold, flip(scale), epsilon);
            theta_temp = theta(j - 1) + dt;
            
            if theta_temp > pi/4
                break
            else
                if j < num_limit
                    theta(j) = theta_temp;
                else
                    error(['The number of the sampled points exceeds the limit of ', ...
                        num2str(num_limit * 4),...
                        '. Please increase the arclength or raise the limit'])
                end
            end
        end
        
        num_pt = j - 1;
        theta = theta(1 : num_pt);
        
        points_fw = angle2points(theta(1 : critical - 1), scale, epsilon);
        points_bw = flip(angle2points(theta(critical : end), flip(scale), epsilon), 2);
        point = [points_fw, [points_bw(2, :); points_bw(1, :)]];
        
        point = [point, flip([-point(1, 1 : num_pt - 1); point(2, 1 : num_pt - 1)], 2), ...
            [-point(1, 2 : end); -point(2, 2 : end)], flip([point(1, 1 : num_pt - 1); ...
            -point(2, 1 : num_pt - 1)], 2)];
        
    end

    function [dt] = dtheta(theta, arclength, threshold, scale, sigma)
        if theta < threshold
            dt = abs((arclength / scale(2) + (theta)^(sigma))^(1 / sigma) ...
                - (theta));
        else
            dt = arclength / sigma * ((cos(theta) ^ 2 * sin(theta) ^ 2) / ...
                (scale(1) ^ 2 * cos(theta) ^ (2 * sigma) * sin(theta) ^ 4 + ...
                scale(2) ^ 2 * sin(theta) ^ (2 * sigma) * cos(theta) ^ 4))^(1 / 2);
        end
    end

    function [point] = angle2points(theta, scale, sigma)
        point = zeros(2, size(theta, 2));
        point(1, :) = scale(1) .* sign(cos(theta)) .* abs(cos(theta)).^sigma;
        point(2, :) = scale(2) .* sign(sin(theta)) .* abs(sin(theta)).^sigma;
    end

end