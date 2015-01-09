function [usr_par] = problemsetup(nxint, nyint, myRHS, myDiffusivity)

%  function [u] = problemsetup(nxint, nyint, nPDEs, generateRHS)
%
%  PURPOSE: Generates problem-specific constant quantities, which require
%           some computational work, but can be computed once and reused,
%           thus speeding up the overall computation. These can include
%           various types of data, and depend entirely on the problem
%           (for example, discretization data).
%           Also, all domain decomposition constructs are computed.
%
%  where y is the solution of the PDE:
%
%  - div(k grad(u))  = f              in Omega
%                  u = u_D            on Gamma_D
%     (k grad(u))'*n = g              on Gamma_N
%
%  The problem domain Omega is the square (xmin,xmax)x(ymin,ymax).
%
%  u     - state
%  k     - material parameter
%  u_exp - observed state (experiments)
%  f     - source term
%
%  INPUT:   nxint          integer
%                          number of subintervals in x-direction
%
%           nyint          integer
%                          number of subintervals in y-direction
%
%           myRHS          functor
%                          right-hand side functor
%
%           myDiffusivity  functor
%                          material property (diffusivity) functor
%
%  OUTPUT:  usr_par     contains all input parameters, as well as additional
%                       computed quantities
%
%  NOTE: Please read in-code comments for a precise explanation of the
%        computed data fields of the usr_par structure.
%
%  AUTHORS: Miguel Aguilo
%           Denis Ridzal
%           (Sandia National Labs)
%           March 30, 2011

%%% General input data
spaceDim  = 2;           % physical spatial dimensions
cellType  = 'Triangle';  % cell type
nVert     = 3;           % number of cell vertices
cubDegree = 3;           % max. degree of the polynomial that can be represented by the basis
numFields = 3;           % number of fields (i.e. number of basis functions)
xmin = -1; xmax = 1;     % min. and max. dimensions along the x-coordinate
ymin = -1; ymax = 1;     % min. and max. dimensions along the y-coordinate

usr_par.spaceDim      = spaceDim;
usr_par.cellType      = cellType;
usr_par.nVert         = nVert;
usr_par.cubDegree     = cubDegree;
usr_par.numFields     = numFields;

%%%%%%%%%% Generate computational mesh on [xmin,xmax]x[ymin,ymax]
[mesh] = RectGrid( xmin, xmax, ymin, ymax, nxint, nyint, cellType);

% Store the number of vertices and cells in the grid
nVertGrid = max( max(mesh.t(:,1:3)) );
numCells  = size(mesh.t, 1);
usr_par.mesh	  = mesh;
usr_par.numCells  = numCells;
usr_par.nVertGrid = nVertGrid;

% Generate degree-of-freedom maps in Intrepid-compatible format
[cellNodes, iIdx, jIdx, iVecIdx] = dofmaps(mesh, nVert, numCells);

% Store DoF data.
usr_par.cellNodes = cellNodes;
usr_par.iIdx      = iIdx;
usr_par.jIdx      = jIdx;
usr_par.iVecIdx   = iVecIdx;

%%%%%%%%%%% Generate reusable finite element data.
[usr_par] = femdata(usr_par);

%%%%%%%%%%% Generate the right hand side.
f = myRHS(usr_par);
usr_par.f = f;

% Generate the diffusivity (material parameter)
k = myDiffusivity(usr_par);
usr_par.k = k;

end
