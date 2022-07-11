function x0 = superquadricInit(point, T)

t0 = sum(point .* T, 2) / sum(T);
[eigenVector, eigenValue] = EigenAnalysis(point - t0, T);
eigenVector = [eigenVector(:, 1), eigenVector(:, 3), cross(eigenVector(:, 1), eigenVector(:, 3))];
s = 1/2 * sum(eigenValue,'all');
a = sqrt(abs(s - eigenValue(1,1)) *5 / sum(T));
b = sqrt(abs(s - eigenValue(3,3)) *5 / sum(T));
c = sqrt(abs(s - eigenValue(2,2)) *5 / sum(T));
euler0 = rotm2eul(eigenVector);
x0 = [1, 1, [a b c]/4, euler0, t0' 0 0];

end

function [EigenVector, EigenValue] = EigenAnalysis(point, T)
    
A = (point.*sqrt(T)) * (point.*sqrt(T))';
B = sum(diag((point.*sqrt(T))' * (point.*sqrt(T))));
MOI = B * eye(3) - A;
[EigenVector, EigenValue] = eig(MOI);
% EigenVector = flip(EigenVector, 2);


end