function [u] = pdesolve(usr_par)
%
%   pdesolve(usr_par)
%
%   PURPOSE: Solve the following PDE using FEM:
%       
%       - div(k grad(u))  = f              in Omega
%                       u = u_D            on Gamma_D
%          (k grad(u))'*n = g              on Gamma_N
%
%   The problem domain Omega is the square (xmin,xmax)x(ymin,ymax).
%
%   Input:
%           usr_par    contains all input parameters
%  
%   Output:
%   
%           u          FEM solution
%
%   AUTHOR:  Miguel Aguilo
%            Denis Ridzal
%            Sandia National Laboratories
%            March 30, 2011

spaceDim      = usr_par.spaceDim;
nVertGrid     = usr_par.nVertGrid;
numFields     = usr_par.numFields;
numCubPoints  = usr_par.numCubPoints;
numCells      = usr_par.numCells;

%%%%%%%%%%% Initialization of lhs (i.e. solution) vector.
u = zeros(nVertGrid,1);

% grab diffusivity and source term
k = usr_par.k;
f = usr_par.f;

% evaluate material parameter k (diffusion coefficients) at the cubature points
k = k( usr_par.mesh.t');
k_at_cub_points = zeros(numCubPoints, numCells);
intrepid_evaluate(k_at_cub_points, k, ...
    usr_par.transformed_val_at_cub_points);

%%%%%%%%%%% combine transformed gradients with diffusion parameters
k_times_transformed_grad_at_cub_points = zeros(spaceDim, ...
    numCubPoints, numFields, numCells);
intrepid_scalarMultiplyDataField(k_times_transformed_grad_at_cub_points, ...
    k_at_cub_points, usr_par.transformed_grad_at_cub_points);

%%%%%%%%%%% integrate stiffnes matrix
cell_stiffness_matrices = zeros(numFields, numFields, numCells);
intrepid_integrate(cell_stiffness_matrices, ...
    k_times_transformed_grad_at_cub_points, ...
    usr_par.weighted_transformed_grad_at_cub_points, 'COMP_BLAS');

%%%%%%%%%%% grab the integrated rhs 
cell_rhs = f;

%%%%%%%%%%% build global stiffness matrix
cell_stiffness_matrices = reshape(cell_stiffness_matrices, 1, ...
    numel(cell_stiffness_matrices));
stiff_mat = sparse(usr_par.iIdx, usr_par.jIdx, cell_stiffness_matrices);

%%%%%%%%%%% build global rhs vector
cell_rhs = reshape(cell_rhs, 1, numel(cell_rhs));
rhs_mat = sparse(usr_par.iVecIdx, usr_par.iVecIdx, cell_rhs);
b = spdiags(rhs_mat,0);

%%%%%%%%%%% Apply Dirichlet conditions to rhs vector. 
if( ~isempty(u) )
    u(unique(usr_par.dirichlet)) = usr_par.u_dirichlet( ...
        unique(usr_par.dirichlet) );
    b = b - stiff_mat * u;
end

%%%%%%%%%%% Computation of the solution
%u(usr_par.FreeNodes) = stiff_mat(usr_par.FreeNodes,usr_par.FreeNodes) ...
%    \ b(usr_par.FreeNodes);

%%%%%%%%%%% set up and run ML
ml('cleanup');
ml_setup = { ...
  'coarse: type',              'Amesos-KLU', ...
  'ML output',                  0, ...
  'coarse: max size',           32, ...
  'smoother: type',            'symmetric Gauss-Seidel' ...
};
ml_apply = { ...
  'krylov: type',              'fixed point', ...
  'krylov: max iterations',     20, ...
  'krylov: tolerance',          1e-100, ...
  'krylov: output level',	1 ...
};
rr = symrcm(stiff_mat(usr_par.FreeNodes,usr_par.FreeNodes));
[h,oc] = ml('setup', stiff_mat(usr_par.FreeNodes(rr),usr_par.FreeNodes(rr)), ml_setup{:});
u(usr_par.FreeNodes(rr)) = ml(h, stiff_mat(usr_par.FreeNodes(rr),usr_par.FreeNodes(rr)), b(usr_par.FreeNodes(rr)), ml_apply{:});

end
