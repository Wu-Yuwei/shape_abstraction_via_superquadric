function [dist] = distance(X, para)
    % transform pose parameters into R matrix and t vector
    R = eul2rotm(para(6 : 8));
    t = para(9 : 11);
    % align the point cloud to the superquadrics coordinate
    X_c = R' * X - R' * t';
    xc = X_c(1,:);
    yc = X_c(2,:);
    zc = X_c(3,:);
    
    kx = para(12);
    ky = para(13);

    xc = xc ./ (kx/para(3) * zc + 1);
    yc = yc ./ (ky/para(3) * zc + 1);
    X_c = [xc;yc;zc];
    % calulate the radial distance of each point
    r_0 = vecnorm(X_c);
    dist = r_0 .* abs(((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
        ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
        ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ (-para(1) / 2) - 1);
    
%     dist = abs(((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
%         ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
%         ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ (-para(1) / 2) - 1);
end