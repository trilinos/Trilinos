function varargout=ml(varargin)
% Interface to the ML.  Sets Up/Runs/Queries/Cleans up the smoothed
% aggregation algebraic multigrid solver ML.  In general, options
% are identical to the Teuchos parameters found in the ML
% documentation. 
%
% Call Syntax:
% (1) [h,oc]=ml('setup',A,['parameter',value,...]) - sets up ML
%  Input:
%  A               - Matrix to be solved with [sparse matrix]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  h               - The ML handle for the system
%  [oc]            - Operator complexity of the preconditioner
%
%  Some standard ML options include:
%  PDE equations   - number of PDE equations for each node [int]
%  max levels      - maximum number of hierarchy levels [int]
%  aggregation: threshold - drop tolerance for aggregation [double]
%  smoother: sweeps- number of smoother sweeps [int]
%  smoother: type  - type of smoother to use [string]
%                    Examples: Jacobi, Gauss-Seidel, Hiptmair
%                    symmetric Gauss-Seidel, MLS
%  Some options only function for MLAPI:
%  additional candidates - number of adaptive vectors [int]
%  use default null space - use the default MLAPI null space [int]
%  Mlmex has one unique option of its own:
%  mlmex: interface- tells mlmex whether to use MLAPI or ML_Epetra
%                    Values: 'mlapi' or 'epetra'  [string]
%                    Default: 'epetra'  
%
% (2) [h,oc]=ml('setup',Kedge,Knode,Grad,['parameter',value,...]) -
% sets up ML for Maxwell's equations
%  Input:
%  Kedge           - Matrix to be solved with [sparse matrix]
%  Knode           - Auxillary nodal Laplacian [sparse matrix]
%  Grad            - Discrete gradient
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  h               - The ML handle for the system
%  [oc]            - Operator complexity of the preconditioner
%
%
% (3) ml(h,A,b,['parameter',value,...])
%  Input:
%  h               - The ML handle for the system to solve [int]
%  A               - Matrix to be solved with [sparse matrix]
%  b               - The RHS to solve with [vector]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%
%  Some standard ML options include:
%  krylov: tolerance- Tolerance for the solver [double]
%  krylov: max iterations- Maximum number of krylov iterations [int]
%  krylov: type    - Solver to use [string]
%                    Examples: cg, gmres, fixed point
%  Even if the ML_Epetra interface is used, these options are given
%  in the MLAPI style.
%
% (4) ml('cleanup',[h]) - frees allocated memory
%  Input:
%  h               - The ML handle for the system to clean up.
%                   Calling 'cleanup' with no handle cleans up all
%                   the systems.
%
% (5) ml('status',[h]) - prints out status information
%  Input:
%  h               - The ML handle for the system to print status for.
%                   Calling 'status' with no handle prints status
%                   for all the systems.
%
% (6) agg=ml('aggregate',A,['parameter',value,...]) - Does a single
% level of aggregation on the matrix A.
%  Input:
%  A               - Matrix to be solved with [sparse matrix]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  agg             - Vector indicating which aggregate each node it
%                    put in
%  This mode does not store a handle for the aggregation.  Instead
%  it just returns the aggregates.  Use this option if you want to
%  turn mlmex into a glorified interface for your favorite
%  ML-supported aggregation routine.
%
% by: Chris Siefert <csiefer@sandia.gov>

% Version History
% 04/26/2013 - Adding Maxwell support.
% 10/04/2006 - Finally added the aggregation interface
% 08/30/2006 - Added ML_epetra interface functionality.
% 08/15/2006 - Added operator complexity handling.
% 08/07/2006 - Added status checking functionality.
% 08/03/2006 - Added functionality to allow ML to have multiple systems ready
%              to be solved at once.  This entailed adding the system handle stuff.
% 05/22/2006 - Teuchos-friendly version.
% 05/16/2006 - Initial Version.

if(nargin>=2 && strcmp(varargin{1},'setup'))
  % Setup mode
  [out,oc]=mlmex(0,varargin{2:nargin});
  varargout{1}=out;
  if(nargout==2),varargout{2}=oc;end 
elseif(nargin>=3 && isnumeric(varargin{1}) && issparse(varargin{2}))
  % Solve mode
  if(size(varargin{2},1)~=length(varargin{3})), fprintf('ML: Error size mismatch between A + B\n');out=0;
  else [sol,its]=mlmex(1,varargin{:}); end
  varargout{1}=sol;
  if(nargout==2), varargout{2}=its;end
elseif(nargin>=1 && strcmp(varargin{1},'cleanup'))  
  % Cleanup mode
  varargout{1}=mlmex(2,varargin{2:nargin});
elseif(nargin>=1 && strcmp(varargin{1},'status'))  
  % Status mode
  varargout{1}=mlmex(3,varargin{2:nargin});  
elseif(nargin>=2 && issparse(varargin{2}) && strcmp(varargin{1},'aggregate'))
  % Aggregate mode
  varargout{1}=mlmex(4,varargin{2:nargin});
else
  fprintf('ML: Error! Invalid input, mode unspecified or insufficient parameters\n');
end



