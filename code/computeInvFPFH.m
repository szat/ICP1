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

function [invFPFH] = computeInvFPFH(points,normals,kdtree,nbh)
% Addaptation of Rusu's paper on the FPFH point descriptors for scale invariance.

% Arguments:  
%           points  - a set of points m x 3
%          normals  - a set of corresponding normals n x 3    
%           kdtree  - kdtree build for points
%              nbh  - size of the neighborhood around each point
%
% Returns: 
%           invFPFH - the table of histograms describing each points, n x 15

%First compute SPSFnbh = 18;
idxNbh = knnsearch(kdtree, points, 'k', nbh+1);
%As idxNbh(i,1) = i
idxNbh = idxNbh(:,2:nbh+1);

nb_points = length(points);

featA = zeros(nb_points,nbh); 
featP = zeros(nb_points,nbh);
featT = zeros(nb_points,nbh);
featD = zeros(nb_points,nbh);

for j = 1:nbh
    %Compute Darboux frames
    %u = ni
    %v = (pj?pi) x u
    %w = u x v
 
    dp = points(idxNbh(:,j),:) - points;
    distance = vecnorm(dp,2,2); %still scaled
    
    %dp = dp./vecnorm(dp,2,2);
    u = normals;
    v = cross(dp./distance,u);
    w = cross(u, v);

    nj = normals(idxNbh(:,j),:);
    
    %Stats
    %alpha = v . nj 
    %phi = (u . (pj?pi))/||pj ?pi|| 
    %theta = arctan(w . nj,u . nj)

    alpha = dot(v,nj,2);
    phi = dot(u,dp./distance,2);
    theta = atan2(dot(w,nj,2),dot(u,nj,2));
    
    featA(:,j) = alpha;
    featP(:,j) = phi;
    featT(:,j) = theta;
    featD(:,j) = distance;
end
disp("Finished computing features.")
pause(0.0001)

%Make histograms of 5 bins for each feature alpha, phi, theta
%binA = [-1, -1/64, -1/256, 1/256, 1/64, 1];
%binA = binA*max(abs(min(min(featA))),abs(max(max(featA))));
binA = linspace(min(min(featA)),max(max(featA)),6); 
binP = linspace(-1,1,6);
binT = [-1, -1/16, -1/64, 1/64, 1/16, 1];
binT = pi*binT;
% binP = linspace(-1,1,6);
% binT = linspace(-pi,pi,6);

spfA = zeros(nb_points,5);
spfP = zeros(nb_points,5);
spfT = zeros(nb_points,5);

for i = 1:nb_points
    spfA(i,:) = histcounts(featA(i,:),binA);
    spfP(i,:) = histcounts(featP(i,:),binP);
    spfT(i,:) = histcounts(featT(i,:),binT);
end
disp("Finished computing SPF histograms.")
pause(0.0001)

%Since these are histograms, they make sense
%Second compute FPFH, FPFH(p) = SPF(p)+(1/k)*SUM_i^k( SPF(pk)/||p-pk|| )
fpfhA = spfA;
fpfhP = spfP;
fpfhT = spfT;

%Ver 2
%Normalize the distance with the distance mean
featD = featD./mean(featD,2);

for j = 1:nbh
    fpfhA = fpfhA + spfA(idxNbh(:,j),:)./(nbh*featD(idxNbh(:,j)));
    fpfhP = fpfhP + spfP(idxNbh(:,j),:)./(nbh*featD(idxNbh(:,j)));
    fpfhT = fpfhT + spfT(idxNbh(:,j),:)./(nbh*featD(idxNbh(:,j)));
end

disp("Finished computing FPFH histograms.")
pause(0.0001)

invFPFH = [fpfhA, fpfhP, fpfhT];

end

