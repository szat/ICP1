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

function [trans_final, trans_sequence] = computeOptimization(V_fixed, V_moving, matches)
% Implementation of the optimization loop from the FGR paper. 

% Arguments:  
%          V_fixed  - the set of points against which the alignment is done, n x 3
%         V_moving  - the set of points which we align, m x 3
%          matches  - the set of correpondences between V_fixed and V_moving, k x 2
% Returns: 
%      trans_final  - the final 4 x 4 matrix representing the alignment of V_moving to V_fixed
%   trans_sequence  - the sequence of 4 x 4 matrices representing the successive alignments 
%                     from the optimization loop

%I am not totally sure about how these parameters interact with the point sets
DIV_FACTOR = 1.4;%          // Division factor used for graduated non-convexity
USE_ABSOLUTE_SCALE = 0;%	// Measure distance in absolute scale (1) or in scale relative to the diameter of the model (0)
MAX_CORR_DIST = 0.5;%       // Maximum correspondence distance (also see comment of USE_ABSOLUTE_SCALE)
ITERATION_NUMBER = 64;%		// Maximum number of iteration

trans_sequence = zeros(4,4,ITERATION_NUMBER);
trans = eye(4);

V_copy = V_moving;

par = 100; %diameter

for itr = 1:ITERATION_NUMBER
    %Stopping condition of while loop
    if (mod(itr,4) == 0 && par > MAX_CORR_DIST)
        par = par / DIV_FACTOR;
    end
    
    %Initialize Jr and r
    nvariable = 6;
    JTJ = zeros(nvariable,nvariable);
    JTr = zeros(nvariable,1);
    J = zeros(nvariable,1);

    %r = 0, in cpp
    r2 = 0.0;
    s = ones(length(matches),1); %Not sure what this is yet
    
    for c = 1:length(matches)
        ii = matches(c,1);
        jj = matches(c,2);

        p = V_fixed(ii,:)';
        q = V_copy(jj,:)';
        rpq = p-q;
        
        c2 = c;
        temp = par/(dot(rpq,rpq)+par);
        s(c2) = temp * temp;
        
        J = zeros(size(J));
        J(2) = -q(3);
        J(3) = q(2);
        J(4) = -1;
        r = rpq(1);
        JTJ = JTJ + J * J' * s(c2);
        JTr = JTr + J * r * s(c2);
        r2 = r2 + r * r * s(c2);
        
        J = zeros(size(J));
        J(3) = -q(1);
        J(1) = q(3);
        J(5) = -1;
        r = rpq(2);
        JTJ = JTJ + J * J' * s(c2);
        JTr = JTr + J * r * s(c2);
        r2 = r2 + r * r * s(c2);
        
        J = zeros(size(J));
        J(1) = -q(2);
        J(2) = q(1);
        J(6) = -1;
        r = rpq(3);
        JTJ = JTJ + J * J' * s(c2);
        JTr = JTr + J * r * s(c2);
        r2 = r2 + r * r * s(c2);
        
        r2 = r2 + (par * (1.0 - sqrt(s(c2)))^2);
    end
    %figure
    %histogram(s)
    result = (-JTJ)\JTr; %Equation 8
    
    %delta
    aff_mat = rotz(result(3))*roty(result(2))*rotx(result(1));
    aff_mat(1,4) = result(4);
    aff_mat(2,4) = result(5);
    aff_mat(3,4) = result(6);
    
    trans_sequence(:,:,itr) = aff_mat;
    trans = aff_mat * trans; 

    ROT = trans(1:3,1:3);
    TR = trans(1:3,4);

    V_copy = (ROT*(V_moving') + TR)';
end

trans_final = trans;
