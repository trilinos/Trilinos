% function [MAT,NODES,DN]=laplacianfun(NPTS,[BCS])
% 
% Generates a discretized Laplacian operator in any arbitrarily
% large number of dimensions.  
% Input:
% NPTS   - Vector containing the number of points per dimension of
%          the discretization
% [BCS]  - Boundary conditions for each variable.  0 [default]
%          gives dirchlet (to a point not included in our grid).  1
%          gives periodic bc's in that variable. 
% Output: 
% MAT    - Discretized Laplacian
% NODES  - Location of the nodes
% DN     - Number of dirchlet neighbors per node (does not work
%          with Neuman bc's).
%
% by: Chris Siefert <siefert@cse.uiuc.edu>
% Last Update: 07/05/06 <siefert>
function [MAT,NODES,DN]=laplacianfun(NPTS,varargin)

% Check for BC's
OPTS=nargin-1;
if (OPTS >= 1) BCS=varargin{1};
else BCS=0*NPTS;end
if (OPTS == 2) SCALE=varargin{2};
else SCALE=false;end

% Pull Constants
NDIM=length(NPTS);
N=prod(NPTS);

% Sanity check
if (length(NPTS) ~= length(BCS)) fprintf('Error: Size mismatch between NPTS and BCS\n');return;end

% Compute jumps (normal)
JUMP=cumprod(NPTS);JUMP=[1,JUMP(1:NDIM)];

% Compute jumps (periodic)
JP=JUMP(2:NDIM+1)-JUMP(1:NDIM);

% Diagonal
MAT=2*NDIM*speye(N,N);

% Assembly
for I=1:NDIM,
  VEC=repmat([ones(JUMP(I)*(NPTS(I)-1),1);zeros(JUMP(I),1)],N/JUMP(I+1),1);        
  VEC=VEC(1:N-JUMP(I));
  MAT=MAT-spdiags([zeros(JUMP(I),1);VEC],JUMP(I),N,N) - spdiags([VEC;zeros(JUMP(I),1)],-JUMP(I),N,N);
  if(BCS(I)==1)
    VEC=repmat([ones(JUMP(I),1);zeros(JUMP(I)*(NPTS(I)-1),1)],N/JUMP(I+1),1);            
    VEC=VEC(1:N-JP(I));
    MAT=MAT-spdiags([zeros(JP(I),1);VEC],JP(I),N,N) - spdiags([VEC;zeros(JP(I),1)],-JP(I),N,N); 
  end
end

% Nodal Location
NODES=[1:NPTS(1)]'; 
for I=2:NDIM,  
  SZ=size(NODES,1);
  NODES=[repmat(NODES,NPTS(I),1),reshape(repmat(1:NPTS(I),SZ,1),SZ*NPTS(I),1)];  
end

% Dirchlet Neighbors
DEG=sum(abs(MAT)>0,2);
DN=full(max(DEG)-DEG);
