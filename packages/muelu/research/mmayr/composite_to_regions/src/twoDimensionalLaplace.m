% twoDimensionalLaplace.m
%
% Build FE stiffness matrix for 1D Laplace problem with Dirichlet BC on the left
% and the right.
%
function [ A ] = twoDimensionalLaplace(nNodes)

nx = sqrt(nNodes); ny = nx;

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
     
A = sparse(A);

end
