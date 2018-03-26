function [normal,normalf] = compute_normal(V,F)

% compute_normal - compute the normal of a triangulation
%
%   [normal,normalf] = compute_normal(vertex,face);
%
%   normal(i,:) is the normal at vertex i.
%   normalf(j,:) is the normal at face j.
%
%   Copyright (c) 2004 Gabriel Peyré

[V,F] = check_face_vertex(V,F);
V = V';
F = F';

nface = size(F,1);
nvert = size(V,1);
normal = zeros(nvert,3);

% unit normals to the faces, %changed crossp to cross
normalf = crossp( V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:) );
d = sqrt( sum(normalf.^2,2) ); d(d<eps)=1;
normalf = normalf ./ repmat( d, 1,3 );

% unit normal to the vertex
normal = zeros(nvert,3);
for i=1:nface
    f = F(i,:);
    for j=1:3
        normal( f(j),: ) = normal( f(j),: ) + normalf(i,:);
    end
end
% normalize
d = sqrt( sum(normal.^2,2) ); d(d<eps)=1;
normal = normal ./ repmat( d, 1,3 );

% enforce that the normal are outward
v = V - repmat(mean(V,2), 1,3);
s = sum( v.*normal, 1 );
if sum(s>0)<sum(s<0)
    % flip
    normal = -normal;
    normalf = -normalf;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);
