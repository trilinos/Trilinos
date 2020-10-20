% function [sol,resnrm]=jacobi(A,x0,b,omega,nits,[varargin])
%
% Runs SOR-Jacobi iteration on the matrix A to solve
% the equation Ax=b.  Assumes A is sparse.  Uses
% initial guess x0.
% -------------- Parameters ---------------
% A     = Matrix of the system.
% x0    = Initial guess
% b     = RHS of system
% D     = Matrix diagonal
% omega = Over-relaxation parameter
% nits  = Iteration cap
% -------------- Outputs ------------------
% sol   = Solution vector from SOR-Jacobi
% resrnm= Norm of the residual from the solution sol
%

%Use jacobi() as "solve" function for MatlabSmoother
% CMS: First time through
function [sol] = jacobi(A, x0, b, D)
omega = 0.9;
nits = 5;
% CMS: Later interface
%function [sol]=jacobi2(A,x0,b,D,omega,nits)

% Initializations
sol = x0;
n = size(A, 1);

% Compute Initial Residual Norm
resnrm(1) = norm(b - A * x0);

for I = 1:nits,
  %disp(['MATLAB: Running Jacobi iteration #', I]);
  % Next iterate
  sol = sol + omega * ((b - A * sol) ./ D);
end
end
