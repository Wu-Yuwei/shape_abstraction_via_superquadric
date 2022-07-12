clear
clc
close all

%% read and prepare data
% demo includes chair && table && lamp
pc = pcread('example/chair.ply'); 
pcd = downSample(pc,3500,3000);
point = double(pcd.Location'); 
idx = randperm(length(point));
point = point(:,idx);
V = (max(point(1, :)) - min(point(1, :))) * (max(point(2, :)) ...
    - min(point(2, :))) * (max(point(3, :)) - min(point( 3, :)));
scale = (200/V)^(1/3);
point = point * scale;

% show original point cloud
figure()
showPoints(point) 

%% run
tic
[theta,sigma,Z] = superquadricSegment(point);
toc





