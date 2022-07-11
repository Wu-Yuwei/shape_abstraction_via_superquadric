clear
clc
close all

%% read data
pc = pcread('example/table.ply');
point = double(pc.Location'); 
idx = randperm(length(point));
point = point(:,idx);
V = (max(point(1, :)) - min(point(1, :))) * (max(point(2, :)) ...
    - min(point(2, :))) * (max(point(3, :)) - min(point( 3, :)));
scale = (200/V)^(1/3);
point = point * scale;
figure()
showPoints(point) 

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
para.sigma2 = 0;
para.realtime_rendering = false;

%% run
tic
[theta,sigma,Z,cost,point] = superquadricSegment(p0, alpha, K, point, T, mergeThre, para);
toc





