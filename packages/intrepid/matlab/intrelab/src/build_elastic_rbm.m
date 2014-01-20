%function [RBM, NDOF]=build_elastic_rbm(NODES)
% Given a set of nodes in 1, 2 or 3 dimensions construct the
% rigid-body modes for the (elastic) structure
% by: Chris Siefert
% Modified by Denis Ridzal on 10/24/2012:
%   - fixed performance issue with sparse vs. zeros
%   - returned dimension and number of rigid body modes
%
function [RBM, dim, numRBM]=build_elastic_rbm(NODES)
% Constants
[N,dim]=size(NODES);
NDOF=[1,3,6];
numRBM=NDOF(dim);

% Allocs
%RBM=sparse(dim*N,numRBM);
RBM=zeros(dim*N,numRBM);

% Translational DOFs / INDICES
for I=1:dim,
  IDX{I}=I:dim:dim*N;
  RBM(IDX{I},I)=ones(N,1); 
end

% Recenter nodes
CTR=sum(NODES) / N;
CNODES = NODES - repmat(CTR,N,1);

% Rotational DOF:  Thanks to Farhat, Pierson and Lesoinne 2000,
% equation (52), you know, once I've managed to make sense out of
% the blasted thing.
if(dim>=2)
  % Rotate in X-Y Plane (around Z axis): [-y ;x];
  RBM(IDX{1},dim+1)=-CNODES(:,2);
  RBM(IDX{2},dim+1)= CNODES(:,1);  
end
if(dim==3)
  % Rotate in Y-Z Plane (around X axis): [-z;y]
  RBM(IDX{2},dim+2)=-CNODES(:,3);
  RBM(IDX{3},dim+2)= CNODES(:,2);    

  % Rotate in X-Z Plane (around Y axis): [z ;-x]
  RBM(IDX{1},dim+3)= CNODES(:,3);
  RBM(IDX{3},dim+3)=-CNODES(:,1);      
end
