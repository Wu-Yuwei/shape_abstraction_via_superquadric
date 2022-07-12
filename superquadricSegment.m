function [theta,sigma,Z] = superquadricSegment(point)

% Written by Yuwei Wu @ NUS
%            Weixiao Liu @ JHU, NUS
% -------------------------------------------------------------------------
% DESCRIPTION: This algorithm solves for the shape abstraction of a given
%              point cloud via nonparametric Bayesian model. By modeling
%              the generation of points as observations sampled from an 
%              infinite mixture of Gaussian Superquadric Taper Models (GSTM). 
%              Our approach formulates the abstraction as a clustering 
%              problem, in which: 1) each point is assigned to a cluster via 
%              the Chinese Restaurant Process (CRP); 2) a primitive representation 
%              is optimized for each cluster, and 3) a merging post-process 
%              is incorporated to provide a concise representation.
%
% INPUT: point - point cloud array (3 x N)
%
% OUTPUT: theta - fitted superquadrics parameters
%         sigma - noise factor of each fitted superquadric
%         Z     - cluster assignment of each point      
% -------------------------------------------------------------------------

%% default parameter
p0 = 0.2;
alpha = 0.5;
K = 30;
T = 30;
mergeThre = 5e-3;
para.iterMax = 4;
para.iterMin = 2;
para.tolerance = 5e-3;
para.relative_tolerance = 0.1;
para.iterLSQ_max = 2;
para.max_switch = 2;
para.adaptive_upper = 1;
% change to true if want to visualize each iteration
para.realtime_rendering = false;

%% setup
% average distance
point = unique(point', 'rows');
[~, averageDist] = knnsearch(point,point,'K',6);
averageDist = averageDist(:,2:end);
averageDist = mean(averageDist,'all');
point = point';

% init
N = size(point,2);
Z = zeros(N,K);
fixedZ = zeros(N,K);
threshold = 20;

%% initilization
sigma = zeros(1,K);
theta = zeros(K,13);
cost = inf(1,K);
Ik = eye(K);

% Firstly segment the point cloud into K clusters by K-means
rng(1)
pIdx = kmeans(point',K, 'Start','uniform', 'Distance','cityblock');
for i = 1:K
    idx = pIdx == i;
    Z(idx,:) = repmat(Ik(i,:),sum(idx),1);
    sigma(i) = rand;
    p = point(:,idx');
    theta(i,:) = superquadricInit(p,ones(1,size(p,2)));
end

%% optimization-based Gibbs sampling
for iter = 1:T
    
    %% separate && split
    [cost, theta, sigma, Z, fixedZ, K, n, Ik] = ...
        split_separate(point, cost, theta, sigma, Z, fixedZ, iter, K, averageDist, threshold, Ik, para);

    %% sample Z
    flag = 0;
    [cor,~] = correspondence_new(point, theta, sigma.^2);
    P = point;
    while ~isempty(P)
        n = sum(Z);
        p = P(:,1);
        pp = find(sum(abs(point-p))==0);
        assert(length(pp)==1)
        nj = Z(pp,:);
        nj_fixed = fixedZ(pp,:);
        if sum(nj_fixed) > 0
            assert(sum(nj_fixed)==1)
            Z(pp,:) = fixedZ(pp,:);
            P(:,1) = [];
            continue
        end
        if flag == 0
            weight = [(n-nj)/(alpha+N-1).*cor(pp,:) alpha/(alpha+N-1)*p0];
        else
            weight = [(n(1:end-1)-nj(1:end-1))/(alpha+N-1).*cor(pp,:) (alpha+n(end))/(alpha+N-1)*p0];
        end
        z = randsample(K+1,1,true,weight);
        
        P(:,1) = [];
    
        if z == K+1 && flag == 0
            Z = [Z zeros(N,1)];
            fixedZ = [fixedZ zeros(N,1)];
            Ik = eye(K+1);
            flag = 1;
            theta = [theta;zeros(1,13)];
            sigma = [sigma rand];
            cost = [cost inf];
        end
        Z(pp,:) = Ik(z,:);
    
        if flag == 1
            newPoint = point(:,(Z(:,K+1)==1)');
            if length(newPoint) > threshold
                [labels,~] = pcsegdist(pointCloud(newPoint'), averageDist*2);
                noCluster = mode(labels);
                cIdx = labels == noCluster;
                if sum(cIdx) > threshold
                    K = K + 1;
                    p = newPoint(:,cIdx');
                    theta(K,:) = superquadricInit(p,ones(1,size(p,2)));
                    flag = 0;
                    [cor_new,~] = correspondence_new(point, theta(K,:), sigma(K)^2);
                    cor = [cor cor_new];
                    np = newPoint(:,~cIdx');
                    P = [P, np];
                end
            end
        end
    
    end
    
    n = sum(Z);
    assert(sum(n)==N)
    K = size(Z,2);
    
    if i <= 10
        A = [];
        for j = length(n):-1:1
            if n(j) < 13
                theta(j,:) = [];
                sigma(j) = [];
                cost(j) = [];
                Z(:,j) = [];
                fixedZ(:,j) = [];
                K = K - 1;
                A = [A,j];
            end
        end
        for j = A
            n(j) = [];
        end
        K = size(Z,2);
        Ik = eye(K);
        
    else
        P = [];
        A = [];
        NIdx = [];
        for j = 1:K
            clusterPoint = point(:,(Z(:,j)==1)');
            [labels,~] = pcsegdist(pointCloud(clusterPoint'), averageDist*2);
            noCluster = mode(labels);
            cIdx = labels == noCluster;
            if sum(cIdx) <= threshold
                P = [P,clusterPoint];
                A = [A,j];
            else
                P = [P,clusterPoint(:,~cIdx)];
            end
        end
        
        for j = flip(A)
            n(j) = [];
            theta(j,:) = [];
            sigma(j) = [];
            cost(j) = [];
            Z(:,j) = [];
            fixedZ(:,j) = [];
            K = K - 1;
        end
        Ik = eye(K);
        
        for j = 1:size(P,2)
            [cor,~] = correspondence_new(P, theta, sigma.^2);
            n = sum(Z);
            p = P(:,j);
            pp = find(sum(abs(point-p))==0);
            nj = Z(j,:);
            nj_fixed = fixedZ(pp,:);
            if sum(nj_fixed) > 0
                assert(sum(nj_fixed)==1)
                Z(pp,:) = fixedZ(pp,:);
                continue
            end
            weight = (n-nj)/(alpha+N-1).*cor(j,:);
            if sum(weight) == 0
                continue
            end
            z = randsample(K,1,true,weight);
            Z(pp,:) = Ik(z,:);
        end
        
    end

    %% optimize parameters of each superquadric
    for j = 1:K
        if iter > 15
            idx = Z(:,j) == 1;
            p = point(:,idx');
            [labels,~] = pcsegdist(pointCloud(p'), averageDist*2);
            noCluster = mode(labels);
            cIdx = labels == noCluster;
            p = p(:,cIdx');
            x0 = superquadricInit(p, ones(1,sum(cIdx)));
            para.iterMax = floor(iter/2);
            [x, D, ~] = superquadricFitting(p, para, x0);
            theta(j,:) = x;
            cost(1,j) = D;
            notp = find(idx ==1);
            Z(notp(~cIdx),j) = 0;
            fixedZ(notp(~cIdx),j) = 0;
        else
            idx = Z(:,j) == 1;
            x0 = superquadricInit(point(:,idx'), ones(1,sum(idx)));
            [x, D, ~] = superquadricFitting(point(:,idx'), para, x0);
            theta(j,:) = x;
            cost(1,j) = D;
        end
    end
    
    c = cost ./ sum(Z);
    for fixedIdx = find(c<=5e-3)
        fixedZ(:,fixedIdx) = Z(:,fixedIdx);
    end    
    
    %% sample sigma
    for j = 1:K
        a = (n(j)-1)/2;
        b = 2/cost(j);
        r = gamrnd(a,b);
        sigma(1,j) = sqrt(1/r);
    end        
    
    %% merge
    if mod(iter,15) == 0
        [n,Z,fixedZ,theta,sigma,K,Ik,cost] = superquadricMerge(Z,iter,fixedZ,K,theta,sigma,cost,point,para,mergeThre);
    end
    
    if para.realtime_rendering
        visualization(Z, K, theta, point, 2, 3, [0 90])
    end
    
    

end
[n,Z,fixedZ,theta,sigma,K,Ik,cost] = superquadricMerge(Z,iter,fixedZ,K,theta,sigma,cost,point,para, mergeThre);

% visualize result
visualization(Z, K, theta, point, 2, 3, [0 90])



