function solnorm = m2ml_test

addpath ../../mesh

fprintf('\nTesting Intrepid and ML functionality ...\n');

spaceDim  = 2;
cellType  = 'Triangle';
nVert     = 3;
cubDegree = 2;

xmin=-1; xmax=2;
ymin=-2; ymax=4;
%nx=5;
%ny=5;
nx=300;
ny=500;

%%%%%%%%%%% get mesh data
fprintf('%50s', 'get mesh data:  ');
tic
[mesh] = RectGrid(xmin, xmax, ymin, ymax, nx, ny, cellType);
toc
pause(0.01)
fprintf('%50s', 'build index arrays:  ');
tic
iIdx_tmp = zeros(3*size(mesh.t,1),size(mesh.t,2));
jIdx_tmp = zeros(size(mesh.t,1),3*size(mesh.t,2));
for i=1:nVert
  iIdx_tmp(i:nVert:numel(mesh.t),:) = mesh.t;
  for j=1:nVert
    jIdx_tmp(:, (i-1)*nVert+j) = mesh.t(:,i);
  end
end
iIdx = reshape(iIdx_tmp', 1, numel(iIdx_tmp));
jIdx = reshape(jIdx_tmp', 1, numel(jIdx_tmp), 1);
clear iIdx_tmp;
clear jIdx_tmp;
toc
pause(0.01)

cellNodesAll = mesh.p(mesh.t',:)';
numCells     = size(mesh.t, 1);
cellNodes    = reshape(cellNodesAll, 2, nVert, numCells); % format to match Intrepid
clear cellNodesAll;

%%%%%%%%%%% evaluate cubature points and weights
numCubPoints = intrepid_getNumCubaturePoints(cellType,cubDegree);
cubPoints    = zeros(spaceDim,numCubPoints);
cubWeights   = zeros(1,numCubPoints);
intrepid_getCubature(cubPoints,cubWeights,cellType,cubDegree);

%%%%%%%%%%% evaluate cell Jacobians
fprintf('%50s', 'evaluate cell Jacobians:  ');
cellJacobians  = zeros(spaceDim,spaceDim,numCubPoints,numCells);
tic
intrepid_setJacobian(cellJacobians,cubPoints,cellNodes,cellType);
toc
pause(0.01)

%%%%%%%%%%% evaluate inverses of cell Jacobians
fprintf('%50s', 'evaluate inverses of cell Jacobians:  ');
cellJacobianInvs = zeros(spaceDim,spaceDim,numCubPoints,numCells);
tic
intrepid_setJacobianInv(cellJacobianInvs,cellJacobians);
toc
pause(0.01)

%%%%%%%%%%% evaluate determinants of cell Jacobians
fprintf('%50s', 'evaluate determinants of cell Jacobians:  ');
cellJacobianDets  = zeros(numCubPoints,numCells);
tic
intrepid_setJacobianDet(cellJacobianDets,cellJacobians);
toc
pause(0.01)

%%%%%%%%%%% evaluate basis (value, gradient)
numFields = 3;
val_at_cub_points = zeros(numCubPoints,numFields);
grad_at_cub_points = zeros(spaceDim,numCubPoints,numFields);
intrepid_getBasisValues(val_at_cub_points, cubPoints, 'OPERATOR_VALUE', cellType, 1);
intrepid_getBasisValues(grad_at_cub_points, cubPoints, 'OPERATOR_GRAD', cellType, 1);

%%%%%%%%%%% compute cell measures
fprintf('%50s', 'compute cell measures:  ');
weighted_measure = zeros(numCubPoints,numCells);
tic
intrepid_computeCellMeasure(weighted_measure, cellJacobianDets, cubWeights);
toc
pause(0.01)

%%%%%%%%%%% transform gradients
fprintf('%50s', 'transform gradients:  ');
transformed_grad_at_cub_points = zeros(spaceDim,numCubPoints,numFields,numCells);
tic
intrepid_HGRADtransformGRAD(transformed_grad_at_cub_points, cellJacobianInvs, grad_at_cub_points);
toc
pause(0.01)

%%%%%%%%%%% transform values
fprintf('%50s', 'transform values:  ');
transformed_val_at_cub_points = zeros(numCubPoints,numFields,numCells);
tic
intrepid_HGRADtransformVALUE(transformed_val_at_cub_points, val_at_cub_points);
toc
pause(0.01)

%%%%%%%%%%% combine transformed gradients with measures
fprintf('%50s', 'combine transformed gradients with measures:  ');
weighted_transformed_grad_at_cub_points = zeros(spaceDim,numCubPoints,numFields,numCells);
tic
intrepid_multiplyMeasure(weighted_transformed_grad_at_cub_points, weighted_measure, transformed_grad_at_cub_points);
toc
pause(0.01)

%%%%%%%%%%% integrate
fprintf('%50s', 'integrate:  ');
cell_stiffness_matrices = zeros(numFields,numFields,numCells);
tic
intrepid_integrate(cell_stiffness_matrices, ...
                   transformed_grad_at_cub_points, ...
                   weighted_transformed_grad_at_cub_points, ...
                   'COMP_BLAS');
toc
pause(0.01)

%%%%%%%%%%% build global stiffness matrix
fprintf('%50s', 'build global stiffness matrix:  ');
tic
cell_stiffness_matrices = reshape(cell_stiffness_matrices, 1, numel(cell_stiffness_matrices));
stiff_mat = sparse(iIdx, jIdx, cell_stiffness_matrices);
toc
pause(0.01)


%%%%%%%%%%% set up and run ML
ml('cleanup');
ml_setup = { ...
  'coarse: type',              'Amesos-KLU', ...
  'ML output',                  1, ...
  'coarse: max size',           128, ...
  'smoother: type',            'symmetric Gauss-Seidel' ...
};
ml_apply = { ...
  'krylov: type',              'fixed point', ...
  'krylov: max iterations',     20, ...
  'krylov: tolerance',          1e-100, ...
  'krylov: output level',	1 ...
};
stiff_mat = stiff_mat + 1e-6*speye(size(stiff_mat));
rr = symrcm(stiff_mat);
[h,oc] = ml('setup', stiff_mat(rr,rr), ml_setup{:});
x = (1/nx/ny)*ones(size(stiff_mat,1), 1);
tic
y(rr) = ml(h, stiff_mat(rr,rr), x(rr), ml_apply{:});
toc

solnorm = norm(y);

fprintf('Done testing Intrepid and ML functionality.\n\n');

end % m2ml_test
