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

function [matches2] = computeMatches(hist1,hist2)
% Implementation of the Reciprocity Test and the Tuple Test from the FGR paper. 

% Arguments:  
%            hist1  - a set of of description histograms for points, n x 15 usually
%            hist2  - a set of of description histograms for points, m x 15 usually
%
% Returns: 
%          matches  - a set of matches between the points represented in hist1 and hits2, k x 2

%Correspondances
kdtree_FPFHinv_init = KDTreeSearcher(hist1);
idx_init = knnsearch(kdtree_FPFHinv_init, hist2, 'k', 3); %returns indices in FPFH1

kdtree_FPFHinv_tform = KDTreeSearcher(hist2);
idx_tform = knnsearch(kdtree_FPFHinv_tform, hist1, 'k', 3); %returns indices in FPFH1

%Reciprocity: idx_init has the size of tform since idx_init are searches in init 
%tform(i) corresponds to init(idx_init(i))
matches = zeros(length(hist2),2);
linear = (1:length(hist2))';
logic = idx_tform(idx_init(linear,1),1) == linear;
matches(logic,2) = linear(logic);
matches(logic,1) = idx_init(logic,1);

%for i = 1:length(invFPFH2)
%   if(idx_tform(idx_init(i,1),1) == i)
%       matches(i,2) = i;
%       matches(i,1) = idx_init(i,1);
%   end
%end

matches(~any(matches,2),:) = [];  %rows
%Now V_init(matches(i,1),:) corresponds to V_tform(matches(i,1),:)

%Replace Tuple Test by Lowe's Ratio Test
%pts_init = V_init(matches(:,1),:);
%pts1_tform = V_tform(idx_tform(matches(:,1),1),:); %are the first picks 
%pts2_tform = V_tform(idx_tform(matches(:,1),2),:); %are the second picks

fts_init = hist1(matches(:,1),:);
fts1_tform = hist2(idx_tform(matches(:,1),1),:);
fts2_tform = hist2(idx_tform(matches(:,1),2),:);

ratio = 0.9;
logic_fts = vecnorm(fts_init - fts1_tform, 2,2) < ratio * vecnorm(fts_init - fts2_tform, 2,2); 
matches = matches(logic_fts,:);

if(length(matches) > 10000)
   matches = datasample(matches, 5000,'Replace',false); 
end

matches2 = matches;
end

