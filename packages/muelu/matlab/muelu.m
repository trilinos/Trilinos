function varargout = muelu(varargin)
%
% Call Syntax:
% (1) [h]=ml('setup',A,['parameter',value,...]) - setup
%  Input:
%  A               - Matrix to be solved with [sparse matrix]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  h               - The MueMex handle for the system
%
%  Muemex has a parameter to control APIs used: "Problem Type"
%  Can be "epetra unprec", "epetra", "tpetra", "tpetra double" or "tpetra complex" 
%	Note: Given a "Problem Type" of "tpetra", Muemex will infer real or complex from the matrix that is passed
%  All but "epetra unprec" use MueLu as a right preconditioner
%  All use Belos PseudoBlockGmres as the solver
%
% (2) [x,its]=ml(h,A,b,['parameter',value,...])
%  Input:
%  h               - The ML handle for the system to solve [int]
%  A               - Matrix to be solved with [sparse matrix]
%  b               - The RHS to solve with [vector]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  x               - Solution vector
%  its             - Number of iterations taken.
%
% (3) [x,its]=ml(h,b,['parameter',value,...])
%  Input:
%  h               - The ML handle for the system to solve [int]
%  b               - The RHS to solve with [vector]
%  parameter       - Teuchos parameter [string]
%  value           - Value for Teuchos parameter [???]
%  Output:
%  x               - Solution vector
%  its             - Number of iterations taken.
%
%  In this case the original 'A' matrix with which this was
%  constructed is now the iteration matrix.
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

if(nargin >= 2 && strcmp(varargin{1},'setup')),
	% Setup mode = 0
	[out, oc] = muemex(0, varargin{2:nargin});
	varargout{1} = out;
	if(nargout == 2),
		varargout{2} = oc;
	end
elseif(nargin >= 3 && isnumeric(varargin{1}) && issparse(varargin{2})),
	% Solve mode (newmatrix)
	if(size(varargin{2}, 1) ~= length(varargin{3})), 
		fprintf('ML: Error size mismatch between A + B\n');
		out = 0;
	else
		[sol, its] = muemex(6, varargin{:});
	end
	varargout{1}=sol;
	if(nargout == 2),
		varargout{2} = its;
  	end
elseif(nargin >= 2 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && (nargin == 2 || ~isnumeric(varargin{3}))),
	% Solve mode (reuse) = 1
	[sol, its] = muemex(1, varargin{:});
	varargout{1} = sol;
	if(nargout == 2),
		varargout{2} = its;
	end
	elseif(nargin >= 1 && strcmp(varargin{1}, 'cleanup')),
	% Cleanup mode = 2
	varargout{1} = muemex(2, varargin{2:nargin});
elseif(nargin>=1 && strcmp(varargin{1}, 'status')),
	% Status mode = 3
	varargout{1} = muemex(3, varargin{2:nargin});  
elseif(nargin>=2 && issparse(varargin{2}) && strcmp(varargin{1}, 'aggregate')),
	% Aggregate mode = 4
	varargout{1} = muemex(4, varargin{2:nargin});
else
	fprintf('Usage:\n');
	fprintf('[h, oc] = muelu(''setup'', A) to setup a problem\n');
	fprintf('muelu(''status''[, probID]) to get status of all problems, or a specific problem.\n');
	fprintf('[x, its] = muelu(h, A, b[, paramName, paramValue, ...]) to solve problem #h\n');
	fprintf('[x, its] = muelu(h, b[, paramName, paramValue, ...]) to solve problem #h with loaded matrix\n');
	fprintf('muelu(''cleanup''[, id]) to free memory associated with all problems, or a specific one.\n');
end
