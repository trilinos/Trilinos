function [usr_par] = femdata(usr_par)
%
%  femdata(usr_par)
%
%  PURPOSE: Generates additional problem-specific FEM quantities, 
%           which require some computational work, but can be computed once 
%           and reused, thus speeding up the overall computation.
%
%  Input:
%           usr_par    contains all input parameters, as well as additional
%                      computed quantities
%  
%  Output:
%   
%           usr_par    contains all previous input parameters, as well as 
%                      additional computed quantities; however, additional
%                      problem-specific constant quantities were added
%
%  AUTHOR:  Miguel Aguilo
%           Denis Ridzal
%           Sandia National Laboratories
%           March 30, 2011

spaceDim      = usr_par.spaceDim;
cellType      = usr_par.cellType;
numFields     = usr_par.numFields;
nVertGrid     = usr_par.nVertGrid;
numCells      = usr_par.numCells;
cubDegree     = usr_par.cubDegree;

%%%%%%%%%%% evaluate cubature points and weights
numCubPoints = intrepid_getNumCubaturePoints(cellType, cubDegree);
cubPoints    = zeros(spaceDim, numCubPoints);
cubWeights   = zeros(1, numCubPoints);
intrepid_getCubature(cubPoints, cubWeights, cellType, cubDegree);
usr_par.numCubPoints = numCubPoints;
usr_par.cubPoints    = cubPoints;
usr_par.cubWeights   = cubWeights;

%%%%%%%%%%%% set boundary markers; mark all boundary edges to be Dirichlet edges
usr_par.mesh.e(:,3) = 1;

%%%%%%%%%%%% Generate boundary data
u_dirichlet = zeros(size(usr_par.mesh.p,1),1);

%%%%%%%%%%% Initialization of free nodes array.
dirichlet = usr_par.mesh.e( (usr_par.mesh.e(:,3) == 1),1:2 );
FreeNodes = setdiff(1:nVertGrid, unique( dirichlet ));

%%%%%%%%%%% evaluate cell Jacobians
cellJacobians  = zeros(spaceDim, spaceDim, numCubPoints, numCells);
intrepid_setJacobian(cellJacobians, usr_par.cubPoints, ...
    usr_par.cellNodes, cellType);
usr_par.cellJacobians = cellJacobians;

%%%%%%%%%%% evaluate inverses of cell Jacobians
cellJacobianInvs = zeros(spaceDim, spaceDim, numCubPoints, numCells);
intrepid_setJacobianInv(cellJacobianInvs, cellJacobians);
usr_par.cellJacobianInvs = cellJacobianInvs;

%%%%%%%%%%% evaluate determinants of cell Jacobians
cellJacobianDets  = zeros(numCubPoints, numCells);
intrepid_setJacobianDet(cellJacobianDets, cellJacobians);
usr_par.cellJacobianDets = cellJacobianDets;

%%%%%%%%%%% evaluate basis (value, gradient)
val_at_cub_points = zeros(numCubPoints, numFields);
grad_at_cub_points = zeros(spaceDim, numCubPoints, numFields);
intrepid_getBasisValues(val_at_cub_points, usr_par.cubPoints, ...
    'OPERATOR_VALUE', cellType, 1);
intrepid_getBasisValues(grad_at_cub_points, usr_par.cubPoints, ...
    'OPERATOR_GRAD', cellType, 1);
usr_par.val_at_cub_points = val_at_cub_points;
usr_par.grad_at_cub_points = grad_at_cub_points;

%%%%%%%%%%% compute cell measures
weighted_measure = zeros(numCubPoints, numCells);
intrepid_computeCellMeasure(weighted_measure, ...
    cellJacobianDets, usr_par.cubWeights);
usr_par.weighted_measure = weighted_measure;

%%%%%%%%%%% transform gradients
transformed_grad_at_cub_points = zeros(spaceDim, numCubPoints, ...
    numFields, numCells);
intrepid_HGRADtransformGRAD(transformed_grad_at_cub_points, ...
    cellJacobianInvs, grad_at_cub_points);
usr_par.transformed_grad_at_cub_points = transformed_grad_at_cub_points;

%%%%%%%%%%% transform values
transformed_val_at_cub_points = zeros(numCubPoints, numFields, numCells);
intrepid_HGRADtransformVALUE(transformed_val_at_cub_points, ...
    val_at_cub_points);
usr_par.transformed_val_at_cub_points = transformed_val_at_cub_points;

%%%%%%%%%%% combine transformed gradients with measures
weighted_transformed_grad_at_cub_points = zeros(spaceDim, ...
    numCubPoints, numFields, numCells);
intrepid_multiplyMeasure(weighted_transformed_grad_at_cub_points, ...
    weighted_measure, transformed_grad_at_cub_points);
usr_par.weighted_transformed_grad_at_cub_points = ...
    weighted_transformed_grad_at_cub_points;

%%%%%%%%%%% combine transformed values with measures
weighted_transformed_val_at_cub_points = zeros(numCubPoints, ...
    numFields, numCells);
intrepid_multiplyMeasure(weighted_transformed_val_at_cub_points, ...
    weighted_measure, transformed_val_at_cub_points);
usr_par.weighted_transformed_val_at_cub_points = ...
    weighted_transformed_val_at_cub_points;

%%%%%%%%%%% integrate stiffness matrix
cell_stiffness_matrices = zeros(numFields, numFields, numCells);
intrepid_integrate(cell_stiffness_matrices, ...
    transformed_grad_at_cub_points, ...
    weighted_transformed_grad_at_cub_points, 'COMP_BLAS');

%%%%%%%%%%% build global stiffness matrix
cell_stiffness_matrices = reshape(cell_stiffness_matrices, 1, ...
    numel(cell_stiffness_matrices));
stiff_mat = sparse(usr_par.iIdx, usr_par.jIdx, cell_stiffness_matrices);

%%%%%%%%%%% integrate mass matrix
cell_mass_matrices = zeros(numFields, numFields, numCells);
intrepid_integrate(cell_mass_matrices, ...
    transformed_val_at_cub_points, ...
    weighted_transformed_val_at_cub_points, 'COMP_BLAS');

%%%%%%%%%%% build global mass matrix
cell_mass_matrices = reshape(cell_mass_matrices, 1, ...
    numel(cell_mass_matrices));
mass_mat = sparse(usr_par.iIdx, usr_par.jIdx, cell_mass_matrices);

%%%%%%%%%%% get computational mesh cubature points physical frame
cubPointsPhysCoord = zeros(spaceDim, numCubPoints, numCells);
intrepid_mapToPhysicalFrame(cubPointsPhysCoord, usr_par.cubPoints, ...
    usr_par.cellNodes, cellType);

usr_par.S           = stiff_mat;
usr_par.M           = mass_mat;
usr_par.nu          = size(FreeNodes, 2);
usr_par.nu_dof      = size(mass_mat, 1);
usr_par.nk          = usr_par.nu_dof;
usr_par.dirichlet   = dirichlet;
usr_par.u_dirichlet = u_dirichlet;
usr_par.FreeNodes   = FreeNodes;
usr_par.cubPointsPhysCoord = cubPointsPhysCoord;

end
