function [n,Z,fixedZ,theta,sigma,K,Ik,cost] = superquadricMerge(Z,iter,fixedZ,K,theta,sigma,cost,point,para, mergeThre)


if iter <=15
    mThre = 1e-3;
else
    mThre = 1e-2;
end

n = sum(Z);
[n,sortIdx] = sort(n,'descend');
theta = theta(sortIdx,:);
sigma = sigma(sortIdx);
cost = cost(sortIdx);
Z = Z(:,sortIdx);
fixedZ = fixedZ(:,sortIdx);

A = [];
for m1 = 1:K-1
    if ismember(m1, A)
        continue
    end
    x1 = theta(m1,:);
    for m2 = m1+1:K
        if ismember(m2, A)
            continue
        end
        x2 = theta(m2,:);
        if abs(x2(1) - x1(1)) > 2
            continue
        end
        p1 = (eul2rotm(theta(m1,6:8)))'*(point(:,Z(:,m1)==1)-theta(m1,9:11)');
        p2 = (eul2rotm(theta(m1,6:8)))'*(point(:,Z(:,m2)==1)-theta(m1,9:11)');
        p3 = [p1,p2];
        V1 = (max(p1(1, :)) - min(p1(1, :))) * (max(p1(2, :)) ...
            - min(p1(2, :))) * (max(p1(3, :)) - min(p1(3, :)));
        V2 = (max(p2(1, :)) - min(p2(1, :))) * (max(p2(2, :)) ...
            - min(p2(2, :))) * (max(p2(3, :)) - min(p2(3, :)));
        V3 = (max(p3(1, :)) - min(p3(1, :))) * (max(p3(2, :)) ...
            - min(p3(2, :))) * (max(p3(3, :)) - min(p3(3, :)));
        if V1 + V2 < 0.7 * V3
            continue
        end

        idx = Z(:,m1) | Z(:,m2) == 1;
        mergedPoint = point(:,idx);
        x0 = superquadricInit(mergedPoint,ones(1,sum(idx)));
        [x, D, ~] = superquadricFitting(mergedPoint, para, x0);
        V1 = 2*x1(1)*x1(2)*x1(3)*x1(4)*x1(5)*beta(x1(1)/2+1,x1(1))*beta(x1(2)/2,x1(2)/2);
        V2 = 2*x2(1)*x2(2)*x2(3)*x2(4)*x2(5)*beta(x2(1)/2+1,x2(1))*beta(x2(2)/2,x2(2)/2);
        V3 = 2*x(1)*x(2)*x(3)*x(4)*x(5)*beta(x(1)/2+1,x(1))*beta(x(2)/2,x(2)/2);
        if V1 + V2 < 0.7 * V3
            continue
        else
            d1 = radial_distance_relative(point(:,Z(:,m1)==1),x1);
            d2 = radial_distance_relative(point(:,Z(:,m2)==1),x2);
            d3 = radial_distance_relative(mergedPoint,x);
            if sum(d3.^2)/sum(idx) <= max(sum(d1.^2)/sum(n(m1)),sum(d2.^2)/sum(n(m2))) || D/sum(idx) < mThre ...
                    || sum(d3.^2)/sum(idx) < mergeThre || D/sum(idx) < max(cost(m1)/n(m1),cost(m2)/n(m2))
%             if D/sum(idx) <= max(cost(m1)/n(m1),cost(m2)/n(m2)) || D/sum(idx) < mThre
                theta(m1,:) = x;
                cost(m1) = D;
                Z(:,m1) = Z(:,m1) | Z(:,m2);
                fixedZ(:,m1) = fixedZ(:,m1) | fixedZ(:,m2);
                sigma(m1) = 1/sqrt(gamrnd((sum(idx)-1)/2, 2/D));
                A = [A, m1,m2];
                n(m1) = sum(idx);
            end
        end

    end
end

A = sort(A(2:2:end),'descend');
for j = A
    n(j) = [];
    Z(:,j) = [];
    fixedZ(:,j) = [];
    theta(j,:) = [];
    sigma(j) = [];
    cost(j) = [];
    K = K - 1;
end
Ik = eye(K);
assert(max(sum(fixedZ,2)) <= 1)

