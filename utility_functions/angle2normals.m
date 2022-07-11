function [normals] = angle2normals(theta, scale, sigma)
normals(1, :) = 1 / scale(1) .* sign(cos(theta)) .* abs(cos(theta)) .^ (2 - sigma);
normals(2, :) = 1 / scale(2) .* sign(sin(theta)) .* abs(sin(theta)) .^ (2 - sigma);
end