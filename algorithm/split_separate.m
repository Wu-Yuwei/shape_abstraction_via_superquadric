function [cost, theta, sigma, Z, fixedZ, K, n, Ik] = ...
    split_separate(point, cost, theta, sigma, Z, fixedZ, i, K, averageDist, threshold, Ik, para)

%% separate
n = sum(Z);
c = cost ./ n;
splitIdx = (i >= 5 & i < 15 & c > 1e-3) | (i >= 15 & c > 1e-2);
if sum(splitIdx)
    A = [];
    for mIdx = find(splitIdx==1)
        splitPoint = point(:,Z(:,mIdx)==1);
        [labels,numClusters] = pcsegdist(pointCloud(splitPoint'), averageDist*2);
        if numClusters == 1
            continue
        end
        noCluster = mode(labels);
        cIdx = labels == noCluster;
        subPoint = splitPoint(:,cIdx);
        if sum(cIdx) > threshold
            x0 = superquadricInit(subPoint, ones(1,sum(cIdx)));
            [theta(mIdx,:),D,~] = superquadricFitting(subPoint, para, x0);
            sigma(mIdx) = sqrt(1/gamrnd((sum(cIdx)-1)/2,2/D));
        else
            A = [A,mIdx];
        end
            
    end
    for j = flip(A)
        n(j) = [];
        Z(:,j) = [];
        fixedZ(:,j) = [];
        theta(j,:) = [];
        sigma(j) = [];
        cost(j) = [];
        K = K - 1;
    end
    Ik = eye(K);        
end

%% split
dThre = 0.7;
if i >= 5
    B = [];
    Ds = [];
    for j = 1:K
        [point_fit] = sphericalProduct_sampling_tapered(theta(j,:), averageDist);
        assert(~isempty(point_fit))
        point_fit = unique(point_fit', 'rows');
        [~, distS2P] = knnsearch(point(:,Z(:,j)==1)',point_fit,'K',1);
%             distS2P = mean(distS2P);
        distS2P = mean(distS2P <= 1.1*averageDist);
        Ds = [Ds, distS2P];
        if distS2P < dThre
            B = [B,i];
            splitPoint = point(:,Z(:,j)==1);
            [labels,numClusters] = pcsegdist(pointCloud(splitPoint'), averageDist*2);
            noCluster = mode(labels);
            cIdx = labels == noCluster;
            splitPoint = splitPoint(:,cIdx);
            rng(1)
            pIdx = kmeans(splitPoint',2, 'Start','uniform', 'Distance','cityblock');
            if sum(pIdx==1) > sum(pIdx==2)
                p = splitPoint(:,pIdx==1);
                if length(p) < threshold
                    continue
                end
                x0 = superquadricInit(p, ones(1,sum(pIdx==1)));
                [theta(j,:),D,~] = superquadricFitting(p, para, x0);
                sigma(j) = sqrt(1/gamrnd((sum(pIdx==1)-1)/2,2/D));
            else
                p = splitPoint(:,pIdx==2);
                if length(p) < threshold
                    continue
                end
                x0 = superquadricInit(p, ones(1,sum(pIdx==2)));
                [theta(j,:),D,~] = superquadricFitting(p, para, x0);
                sigma(j) = sqrt(1/gamrnd((sum(pIdx==2)-1)/2,2/D));
            end
        end
    end

end