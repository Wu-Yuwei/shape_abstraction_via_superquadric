function [T,D] = correspondence_new(point, X, sigma2)
    X = X';
    D = zeros(size(point,2),size(X,2));
    for j = 1:size(X,2)
        D(:,j) = radial_distance(point, X(:,j)');
    end
    
    c = 1 ./ (sqrt(2*pi*sigma2));
    T = c .* exp(-D.^2 ./ (2*sigma2));
   
end