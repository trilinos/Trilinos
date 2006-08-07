function out=ml(varargin)
% Interface to the ML/MLAPI.  Sets Up/Runs/Queries/Cleans up the smoothed
% aggregation algebraic multigrid solver ML.  In general, options
% are identical to the Teuchos parameters found in the ML
% documentation. 
%
% Call Syntax:
% (1) h=ml('setup',A,['parameter',value,...]) - sets up ML
%  Input:
%  A               - Matrix to be solved with [sparse matrix]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  h               - The ML handle for the system
%
%  Some standard MLAPI options include:
%  PDE equations   - number of PDE equations for each node [int]
%  max levels      - maximum number of hierarchy levels [int]
%  additional candidates - number of adaptive vectors [int]
%  use default null space - use the default MLAPI null space [int]
%  aggregation: threshold - drop tolerance for aggregation [double]
%  smoother: sweeps- number of smoother sweeps [int]
%  smoother: type  - type of smoother to use [string]
%                    Examples: Jacobi, Gauss-Seidel, Hiptmair
%                    symmetric Gauss-Seidel, MLS
%
% (2) ml(h,A,b,['parameter',value,...])
%  Input:
%  h               - The ML handle for the system to solve [int]
%  A               - Matrix to be solved with [sparse matrix]
%  b               - The RHS to solve with [vector]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%
%  Some standard MLAPI options include
%  krylov: tolerance- Tolerance for the solver [double]
%  krylov: max iterations- Maximum number of krylov iterations [int]
%  krylov: type    - Solver to use [string]
%                    Examples: cg, gmres, fixed point
%
% (3) ml('cleanup',[h]) - frees allocated memory
%  Input:
%  h               - The ML handle for the system to clean up.
%                   Calling 'cleanup' with no handle cleans up all
%                   the systems.
%
% (4) ml('status',[h]) - prints out status information
%  Input:
%  h               - The ML handle for the system to print status for.
%                   Calling 'status' with no handle prints status
%                   for all the systems.
%
%
% by: Chris Siefert <csiefer@sandia.gov>
% Version History
% 08/07/2006  - Added status checking functionality.
% 08/03/2006  - Added functionality to allow ML to have multiple systems ready
%               to be solved at once.  This entailed adding the system handle stuff.
% 05/22/2006  - Teuchos-friendly version  
% 05/16/2006  - Initial Version 

%

if(nargin>=1 && strcmp(varargin{1},'cleanup'))  
  % Cleanup mode
  out=mlmex(2,varargin{2:nargin});
elseif(nargin>=1 && strcmp(varargin{1},'status'))  
  % Status mode
  out=mlmex(3,varargin{2:nargin});  
elseif(nargin>=2 && strcmp(varargin{1},'setup'))
  % Setup mode
  out=mlmex(0,varargin{2:nargin});
elseif(nargin>=3 && isnumeric(varargin{1}) && issparse(varargin{2}))
  % Solve mode
  if(size(varargin{2},1)~=length(varargin{3})), fprintf('ML: Error size mismatch between A + B\n');out=0;
  else out=mlmex(1,varargin{:}); end
else
  fprintf('ML: Error! Invalid input, mode unspecified or insufficient parameters\n');
end



