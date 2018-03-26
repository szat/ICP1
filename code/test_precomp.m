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

options.sub_type = 'sqrt3';
options.spherical = 0;
options.verb = 0;

[V,F] = readOBJ("data/3d_scan/scan.obj");
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
V = V'; 
F = F';
ptCloud = pointCloud(V);
ptCloud = pcdownsample(ptCloud,'gridAverage',0.5);
normals = compute_normal(V,F);
kdtree = KDTreeSearcher(V);
idx = knnsearch(kdtree, ptCloud.Location, 'k', 1);
normals = normals(idx,:);
V = ptCloud.Location;
kdtree = KDTreeSearcher(V);
features = computeInvFPFH(V,normals,kdtree,33);
pcwrite(ptCloud,"data/3d_scan/scan.ply")
save("data/3d_scan/scan_features",'features')
save("data/3d_scan/scan_normals",'normals')

clear
clc

options.sub_type = 'sqrt3';
options.spherical = 0;
options.verb = 0;

[V,F] = readOBJ("data/chair/chair.obj");
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
V = V';
F = F';
ptCloud = pointCloud(V);
ptCloud = pcdownsample(ptCloud,'gridAverage',0.5);
normals = compute_normal(V,F);
kdtree = KDTreeSearcher(V);
idx = knnsearch(kdtree, ptCloud.Location, 'k', 1);
normals = normals(idx,:);
V = ptCloud.Location;
kdtree = KDTreeSearcher(V);
features = computeInvFPFH(V,normals,kdtree,33);
pcwrite(ptCloud,"data/chair/chair.ply")
save("data/chair/chair_features",'features')
save("data/chair/chair_normals",'normals')

clear
clc

options.sub_type = 'sqrt3';
options.spherical = 0;
options.verb = 0;

[V,F] = readOBJ("data/main/main.obj");
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
V = V';
F = F';
ptCloud = pointCloud(V);
ptCloud = pcdownsample(ptCloud,'gridAverage',0.5);
normals = compute_normal(V,F);
kdtree = KDTreeSearcher(V);
idx = knnsearch(kdtree, ptCloud.Location, 'k', 1);
normals = normals(idx,:);
V = ptCloud.Location;
kdtree = KDTreeSearcher(V);
features = computeInvFPFH(V,normals,kdtree,33);
pcwrite(ptCloud,"data/main/main.ply")
save("data/main/main_features",'features')
save("data/main/main_normals",'normals')

clear
clc

options.sub_type = 'sqrt3';
options.spherical = 0;
options.verb = 0;

[V,F] = readOBJ("data/leg/leg.obj");
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
[V,F] = perform_mesh_subdivision(V', F', 1, options);
V = V';
F = F';
ptCloud = pointCloud(V);
ptCloud = pcdownsample(ptCloud,'gridAverage',0.5);
normals = compute_normal(V,F);
kdtree = KDTreeSearcher(V);
idx = knnsearch(kdtree, ptCloud.Location, 'k', 1);
normals = normals(idx,:);
V = ptCloud.Location;
kdtree = KDTreeSearcher(V);
features = computeInvFPFH(V,normals,kdtree,33);
pcwrite(ptCloud,"data/leg/leg.ply")
save("data/leg/leg_features",'features')
save("data/leg/leg_normals",'normals')
