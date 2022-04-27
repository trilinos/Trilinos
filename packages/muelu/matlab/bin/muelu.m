function varargout = muelu(varargin)
%
% Call Syntax:
% (1) [h,oc] = muelu('setup', A, [Coordinates], [, 'parameter', value,...]) - setup
%  Input:
%  A               - Matrix to be solved with (must be sparse)
%  Coordinates     - User-defined coordinates array (optional, dense array)
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
% (4) x = muelu('apply', h, r)
%  Input:
%  h               - The MueMex handle for the system to solve (int)
%  r               - The residual/RHS vector to solve with (vector or multivector)
%  Output:
%  x               - M*r, where M is the preconditioner operator
%
% (5) muelu('cleanup'[, h]) - frees allocated memory
%  Input:
%  h               - The MueMex handle for the system to clean up.
%                   Calling 'cleanup' with no handle cleans up all
%                   the systems.
%
% (6) muelu('status'[, h]) - prints out status information
%  Input:
%  h               - The MueMex handle for the system to print status for.
%                   Calling 'status' with no handle prints status
%                   for all the systems.
%
% (7) data = muelu('get', h, levelID, dataName[, typeHint])
% - get hierarchy data
% Input:
% h                - MueMex handle for problem with the hierarchy
% levelID          - ID of level to get data from (0, 1, etc.)
% dataName         - Name of data field to fetch ('Nullspace', 'P', etc.)
% typeHint         - (Optional) Type of data expected ('matrix', 'scalar', 'lovector', etc).
%                    If not given, will attempt to infer from dataName.
% Output:
% data	   - The data from the hierarchy.
% 
% by: Brian Kelley <bmkelle@sandia.gov>

if(nargin >= 2 && strcmp(varargin{1},'setup')),
    % Setup mode = 0
    [out, oc] = muemex(0, varargin{2:nargin});
    varargout{1} = out;
    if(nargout == 2),
        varargout{2} = oc;
    end
elseif(nargin >= 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})),
    % Solve mode = 1
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
    % Get mode = 5
    varargout{1} = muemex(5, varargin{2:nargin});
elseif(strcmp(varargin{1}, 'apply')),
    % Apply mode = 8
    varargout{1} = muemex(8, varargin{2:nargin});
else
    fprintf('\nUsage:\n');
    fprintf('[h, oc] = muelu(''setup'', A) to setup a problem\n');
    fprintf('muelu(''status''[, probID]) to get status of all problems, or a specific problem.\n');
    fprintf('[x, its] = muelu(h, A, b[, paramName, paramValue, ...]) to solve problem #h\n');
    fprintf('[x, its] = muelu(h, b[, paramName, paramValue, ...]) to solve problem #h with loaded matrix\n');
    fprintf('x = muelu(''apply'', h, r) to apply preconditioner of problem #h to vector r\n');
    fprintf('muelu(''cleanup''[, id]) to free memory associated with all problems, or a specific one.\n');
    fprintf('muelu(''get'', h, levelID, dataName[, typeHint]) to get data from a level of problem #h.\n');
end
end
