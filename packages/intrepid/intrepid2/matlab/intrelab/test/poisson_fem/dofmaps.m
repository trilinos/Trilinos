% function [cellNodes, iIdx, jIdx, iVecIdx] = dofmaps(mesh, nVert, numCells)
%
%  PURPOSE: Given a mesh in point/triangle/edge format, generate
%           degree-of-freedom maps compatible with Intrepid.
%
%  INPUT:             mesh       structure array with the following fields
%
%                       - mesh.p     Real nn x 2 
%                                    array containing the x- and y-
%                                    coordinates of the nodes
%                                        
%
%                       - mesh.t     Integer nt x 3   
%                                    t(i,1:4) contains the indices of the 
%                                    vertices of triangle i. 
%
%                       - mesh.e     Integer nf x 3
%                                    e(i,1:2) contains the indices of the 
%                                    vertices of edge i.
%                                    edge(i,3) contains the boundary marker 
%                                    of edge i.
%                                    Currently set to one.
%                                    e(i,3) = 1  Dirichlet bdry conds are 
%                                                imposed on edge i
%                                    e(i,3) = 2  Neumann bdry conds are 
%                                                imposed on edge i
%                                    e(i,3) = 3  Robin bdry conds are 
%                                                imposed on edge i
%
%                     nVert      number of cell vertices
%
%                     numCells   number of working cells
%
%  OUTPUT:            cellNodes  Integer 3 x numCells
%                                contains all the indices of the vertices
%                                of cell i.
%
%                     iIdx
%
%                     jIdx
%
%                     iVecIdx
%

function [cellNodes, iIdx, jIdx, iVecIdx] = dofmaps(mesh, nVert, numCells)

iIdx_tmp = zeros(3*size(mesh.t,1),size(mesh.t,2));
jIdx_tmp = zeros(size(mesh.t,1),3*size(mesh.t,2));
for i=1:nVert
  iIdx_tmp(i:nVert:numel(mesh.t),:) = mesh.t;
  for j=1:nVert
    jIdx_tmp(:, (i-1)*nVert+j) = mesh.t(:,i);
  end
end
iIdx = reshape(iIdx_tmp', 1, numel(iIdx_tmp));
jIdx = reshape(jIdx_tmp', 1, numel(jIdx_tmp), 1);
clear iIdx_tmp;
clear jIdx_tmp;

% get rhs array index
iVecIdx = reshape(mesh.t',1,numel(mesh.t));

cellNodesAll = mesh.p(mesh.t',:)';

% format to match Intrepid
cellNodes    = reshape(cellNodesAll, 2, nVert, numCells);
clear cellNodesAll;

end
