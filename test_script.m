clear
clc
close all

addpath('algorithm/')
addpath('utility_functions/')
%% read and prepare data
% demo includes chair && table && lamp
pc = pcread('example/chair.ply'); 
pcd = downSample(pc,3500,3000);
point = double(pcd.Location'); 

% show original point cloud
figure()
showPoints(point) 

%% run
tic
[theta,sigma,Z] = superquadricSegment(point);
toc





