function [point] = angle2points(theta, scale, sigma)
point = zeros(2, size(theta, 2));
point(1, :) = scale(1) .* sign(cos(theta)) .* abs(cos(theta)).^sigma;
point(2, :) = scale(2) .* sign(sin(theta)) .* abs(sin(theta)).^sigma;
end