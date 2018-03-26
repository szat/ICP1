% Copyright (c) Adrian Szatmari
% Author: Adrian Szatmari
% Date: 2018-02-03
% License: MIT, patent permitting
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

clear
close all
clc

pc1 = pcread("data/3d_scan/scan.ply");
normals1 = load("data/3d_scan/scan_normals.mat");
features1 = load("data/3d_scan/scan_features.mat");
V1 = pc1.Location;
features1 = features1.features;

pc2 = pcread("data/chair/chair.ply");
normals2 = load("data/chair/chair_normals.mat");
features2 = load("data/chair/chair_features.mat");
V2 = pc2.Location;
features2 = features2.features;

matches = computeMatches(features1, features2);
[trans_final, trans_sequence] = computeOptimization(V1, V2, matches);

%Viz
V_temp = V2;
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(trans_sequence)
    V_temp = (trans_sequence(1:3,1:3,i)*(V_temp') + trans_sequence(1:3,4,i))';
    
    subplot(1,2,1)
    scatter3(V1(1:4:end,1),V1(1:4:end,2),V1(1:4:end,3));
    hold on
    scatter3(V_temp(1:4:end,1),V_temp(1:4:end,2),V_temp(1:4:end,3));
    view(2)
    title("Front View, iter " +  num2str(i));
    hold off

    subplot(1,2,2)
    scatter3(V1(1:4:end,1),V1(1:4:end,2),V1(1:4:end,3));
    hold on
    scatter3(V_temp(1:4:end,1),V_temp(1:4:end,2),V_temp(1:4:end,3));
    view(-180,-90)
    title("Rear View, iter " +  num2str(i));
    hold off

    pause(0.1);
end

%Save Final Result
% V_temp = trans_final(1:3,1:3)*V2'+trans_final(1:3,4);
% V_temp = V_temp';
% ptCloud = pointCloud(V_temp);
% pcwrite(ptCloud,"data/aligned/aligned.ply");
% [V,F] = readOBJ("data/chair/chair.obj");
% V_temp = trans_final(1:3,1:3)*V'+trans_final(1:3,4);
% V_temp = V_temp';
% writeOBJ("data/aligned/aligned.obj",V_temp,F)
