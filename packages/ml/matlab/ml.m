function out=ml(varargin)
% Interface to the ML.  Sets Up/Runs/Cleans up the smoothed
% aggregation algebraic multigrid solver ML.
% Call Syntax:
% (1) ml('setup',A,[numpdes,cutoff,avecs,agg_thresh,ssweeps]) - sets up ML
%     ml('setup',A,['parameter',value,...])
%  Input:
%   numpdes      - number of PDE equations being solved (default: 1)
%   cutoff       - number of levels to use in the heirarchy (default: 10)
%   avecs        - number of Adaptive SA vectors to use (default: 1)
%   agg_thresh   - aggregation threshold (default: 0)
%                  WARNING: if agg_thresh>0, then we must have avecs=0  
%   ssweeps      - number of smoother sweeps to use (default: 3)
%   
%  The second format allows for arbitrary (Teuchos-style) parameter
%  lists.   These parameters are precisely those supported by the
%  MLAPI and get imported right in.
% 
% Special Options:
% mlmex: numpdes - number of PDEs for SA.
% mlmex: levels  - maximum number of heirarchy levels
%
%
% (2) x=ml(A,b,[iters,tolerance,solver]) - solves w/ ML
%     ml(A,b,['parameter',value,...])
%   Input:
%    iters       - maximum number of solver iterations to allow
%    tolerance   - tolerence for the iterative solver
%    solver      - 1 = ML (fixed point)
%                  2 = CG w/ ML preconditioner
%                  3 = GMRES w/ ML preconditioner
%
%  The second format allows for arbitrary (Teuchos-style) parameter
%  lists.   These parameters are precisely those supported by the
%  MLAPI and get imported right in.
%
% (3) ml('cleanup') - frees allocated memory
%
% by: Chris Siefert <csiefer@sandia.gov>
% Initial Version 05/16/06 
% Teuchos-friendly version 05/22/06 

if(nargin==1 && strcmp(varargin{1},'cleanup'))  
  % Cleanup mode
  mlmex(2);
  out=1;
elseif(strcmp(varargin{1},'setup'))
  % Setup mode
  if(nargin==3) 
    opts=varargin{3};
    if(length(opts)>3 &&isnumeric(opts) && opts(4)>0 && opts(3)>0)
      fprintf('Error: Non-zero aggregation threshold and number of adaptive vectors.  One of these must be zero!');      
      out=1;
      return;
    end
  end
  mlmex(0,varargin{2:nargin});
  out=1;
elseif(issparse(varargin{1}))
  % Solve mode
  out=mlmex(1,varargin{:});
else
  fprintf('ML: Error! Invalid input, mode unspecified\n');
end



