function varargout = muelu(varargin)
%
% Call Syntax:
% (1) [h,oc] = muelu('setup', A[, 'parameter', value,...]) - setup
%  Input:
%  A               - Matrix to be solved with (must be sparse)
%  parameter       - Teuchos parameter key (string)
%  value           - Value for Teuchos parameter
%  Output:
%  h               - The MueMex handle for the system
%  oc              - The operator complexity of MueLu hierarchy
%
%  Muemex has a parameter to control APIs used: "Linear Algebra"
%  Can be "epetra unprec", "epetra", "tpetra"
%	Note: "tpetra" supports real and complex scalar matrices
%  All but "epetra unprec" use MueLu as a right preconditioner
%  All use a Belos solver (PseudoBlockGmres is default)
%
% (2) [x,its] = muelu(h, A, b[, 'parameter', value,...])
%  Input:
%  h               - The MueMex handle for the system to solve (int)
%  A               - Matrix to be solved with
%  b               - The RHS to solve with (vector or multivector)
%  parameter       - Teuchos parameter key (string)
%  value           - Value for Teuchos parameter
%  Output:
%  x               - Solution (vector or multivector)
%  its             - Number of iterations taken by Belos
%
% (3) [x,its] = muelu(h, b[, 'parameter', value,...])
%  Input:
%  h               - The MueMex handle for the system to solve (int)
%  b               - The RHS to solve with (vector or multivector)
%  parameter       - Teuchos parameter key (string)
%  value           - Value for Teuchos parameter
%  Output:
%  x               - Solution (vector or multivector)
%  its             - Number of iterations taken by Belos
%
%  In this case the original 'A' matrix with which this was
%  constructed is now the iteration matrix.
%
% (4) muelu('cleanup'[, h]) - frees allocated memory
%  Input:
%  h               - The MueMex handle for the system to clean up.
%                   Calling 'cleanup' with no handle cleans up all
%                   the systems.
%
% (5) muelu('status'[, h]) - prints out status information
%  Input:
%  h               - The MueMex handle for the system to print status for.
%                   Calling 'status' with no handle prints status
%                   for all the systems.
%
% (6) data = muelu('get', h, levelID, dataName[, typeHint])
% - get hierarchy data
% Input:
% h				   - MueMex handle for problem with the hierarchy
% levelID		   - ID of level to get data from (0, 1, etc.)
% dataName		   - Name of data field to fetch ('Nullspace', 'P', etc.)
% typeHint		   - (Optional) Type of data expected ('matrix', 'scalar',
%													   'lovector', etc).
%						If not given, will attempt to infer from dataName.
% Output:
% data			   - The data from the hierarchy.
%
% by: Brian Kelley <bmkelle@sandia.gov>

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
		fprintf('MueLu: Error size mismatch between A & b\n');
		fprintf('Make sure b has same number of rows as A.\n');
		varargout{1} = 0;
		if(nargout == 2),
			varargout{2} = 0;
		end
	else
		[sol, its] = muemex(6, varargin{:});
		varargout{1} = sol;
		if(nargout == 2),
			varargout{2} = its;
  		end
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
elseif(strcmp(varargin{1}, 'get')),
	% Get mode = 8
	varargout{1} = muemex(8, varargin{2:nargin});
else
	fprintf('\nUsage:\n');
	fprintf('[h, oc] = muelu(''setup'', A) to setup a problem\n');
	fprintf('muelu(''status''[, probID]) to get status of all problems, or a specific problem.\n');
	fprintf('[x, its] = muelu(h, A, b[, paramName, paramValue, ...]) to solve problem #h\n');
	fprintf('[x, its] = muelu(h, b[, paramName, paramValue, ...]) to solve problem #h with loaded matrix\n');
	fprintf('muelu(''cleanup''[, id]) to free memory associated with all problems, or a specific one.\n');
end
