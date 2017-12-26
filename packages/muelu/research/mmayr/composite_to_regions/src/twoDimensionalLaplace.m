% twoDimensionalLaplace.m
%
% Build FE stiffness matrix for 1D Laplace problem with Dirichlet BC on the left
% and the right.
%
function [ A ] = twoDimensionalLaplace(nNodes)

nx = sqrt(nNodes); 
ny = nx;

%% create matrix
A = zeros(nNodes,nNodes);
for j=1:ny, for i=1:nx,
     me = (j-1)*nx + i;

     A(me,me) = 8.;
     if i ~=1, 
        A(me,me-1) = -1;
        if j ~= 1, A(me,me-1-nx) = -1; end;
        if j ~=ny, A(me,me-1+nx) = -1; end;
     end
     if i ~=nx, 
        A(me,me+1) = -1;
        if j ~= 1, A(me,me+1-nx) = -1; end;
        if j ~=ny, A(me,me+1+nx) = -1; end;
     end
     if j ~=  1, A(me,me-nx) = -1; end;
     if j ~= ny, A(me,me+nx) = -1; end;
end;end;

%nzsPerRow = spones(A)*ones(nNodes,1);
%Dirs = find(nzsPerRow == 1);
%A(Dirs,Dirs) = speye(length(Dirs),length(Dirs));

%% impose Dirichlet boundary conditions

% top
inds = 1:nx;
A(inds,:) = zeros(length(inds), size(A,2));
A(inds,inds) = eye(length(inds));

% left
inds = 1:nx:nNodes;
A(inds,:) = zeros(length(inds), size(A,2));
A(inds,inds) = eye(length(inds));

% right
inds = nx:nx:nNodes;
A(inds,:) = zeros(length(inds), size(A,2));
A(inds,inds) = eye(length(inds));

% bottom
inds = nNodes-nx:nNodes;
A(inds,:) = zeros(length(inds), size(A,2));
A(inds,inds) = eye(length(inds));

%% Convert to sparse format    
A = sparse(A);

end
