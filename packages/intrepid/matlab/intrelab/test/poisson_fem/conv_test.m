% Convergence test for Poisson problem.

function rate = conv_test

clear all;

set(0, 'defaultaxesfontsize',12,'defaultaxeslinewidth',0.7,...
    'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.7,...
    'defaulttextfontsize',12);

% mesh utilities directory
addpath ../../mesh/

% choose solution plotting
iplot = 0;

fprintf('\nTesting mesh convergence of Poisson solver ...\n');

% set rhs_fn = @rhs_sine;
rhs_fn=@(usr_par)rhs_sine(usr_par);

% set diff_fn = @diff_const;
diff_fn=@(usr_par)diff_const(usr_par);

usr_par = [];

gridsizes = [4:9];
L2_err = [];

for imesh = gridsizes

  clear usr_par;

  % number of intervals in x and y direction on a rectangular grid
  nxint = 2^imesh;
  nyint = nxint;

  % set up PDE
  fprintf('\nSetting up problem and assembling FEM operators ...\n');
  tic
  [usr_par] = problemsetup(nxint, nyint, rhs_fn, diff_fn);
  toc

  % solve PDE
  fprintf('\nSolving PDE ...\n\n');
  tic
  [u] = pdesolve(usr_par);
  toc

  if iplot
    figure(1)
    trisurf(usr_par.mesh.t, usr_par.mesh.p(:,1), usr_par.mesh.p(:,2), ...
            u, 'facecolor','interp')
    view(0,90);
    axis equal
    axis tight
    shading interp;
    title('PDE Solution')
    xlabel('x')
    ylabel('y')
    pause(1)
  end

  %%% Compute error:

  % set up more accurate numerical integration
  cubDegree = 6;
  numCubPoints = intrepid_getNumCubaturePoints(usr_par.cellType, cubDegree);
  cubPoints    = zeros(usr_par.spaceDim, numCubPoints);
  cubWeights   = zeros(1, numCubPoints);
  intrepid_getCubature(cubPoints, cubWeights, usr_par.cellType, cubDegree);
  val_at_cub_points = zeros(numCubPoints, usr_par.numFields);
  intrepid_getBasisValues(val_at_cub_points, cubPoints, 'OPERATOR_VALUE', usr_par.cellType, 1);
  transformed_val_at_cub_points = zeros(numCubPoints, usr_par.numFields, usr_par.numCells);
  intrepid_HGRADtransformVALUE(transformed_val_at_cub_points, val_at_cub_points);
  
  cellMeasures = zeros(numCubPoints, usr_par.numCells);
  cellJacobians = zeros(usr_par.spaceDim, usr_par.spaceDim, numCubPoints, usr_par.numCells);
  cellJacobianInvs = zeros(usr_par.spaceDim, usr_par.spaceDim, numCubPoints, usr_par.numCells);
  cellJacobianDets = zeros(numCubPoints, usr_par.numCells);

  intrepid_setJacobian(cellJacobians, cubPoints, usr_par.cellNodes, usr_par.cellType);
  intrepid_setJacobianInv(cellJacobianInvs, cellJacobians);
  intrepid_setJacobianDet(cellJacobianDets, cellJacobians);  

  intrepid_computeCellMeasure(cellMeasures, cellJacobianDets, cubWeights);

  % evaluate approximate solution at integration points
  approxSolnCoeffs = reshape(u(usr_par.iVecIdx), usr_par.numFields, usr_par.numCells);  % scatter to cells
  approxSolnValues = zeros(numCubPoints, usr_par.numCells);
  intrepid_evaluate(approxSolnValues, ...
                    approxSolnCoeffs, ...
                    transformed_val_at_cub_points);

  % evaluate true solution at integration points
  % (first map reference integration points to physical frame)
  cubPointsPhysCoords = zeros(usr_par.spaceDim,numCubPoints,usr_par.numCells);
  intrepid_mapToPhysicalFrame(cubPointsPhysCoords, ...
                              cubPoints, ...
                              usr_par.cellNodes, ...
                              usr_par.cellType);
  trueSolnValues = zeros(numCubPoints, usr_par.numCells);
  trueSolnValues = soln_sine(cubPointsPhysCoords);

  % evaluate difference
  diffSolnValues = trueSolnValues-approxSolnValues;

  % evaluate measure-weighted difference
  diffSolnValuesWeighted = zeros(numCubPoints, usr_par.numCells);
  intrepid_scalarMultiplyDataData(diffSolnValuesWeighted, cellMeasures, diffSolnValues);

  % compute L2 error
  err = zeros(1,usr_par.numCells);
  intrepid_integrate(err, diffSolnValues, ...
    diffSolnValuesWeighted, 'COMP_BLAS');
  L2_err = [L2_err sqrt(sum(err))];

end

L2_ratio  = [0 L2_err(1:end-1)]./L2_err;
L2_rate   = 0;
L2_rate_c = 0;

for i=2:length(L2_err)
  p = polyfit(log(1./(2.^gridsizes(i-1:i))), log(L2_err(i-1:i)), 1);
  L2_rate = [L2_rate p(1)];
end

for i=2:length(L2_err)
  p = polyfit(log(1./(2.^gridsizes(1:i))), log(L2_err(1:i)), 1);
  L2_rate_c = [L2_rate_c p(1)];
end

fprintf('\n    mesh         L2 error       ratio          rate_two_grid        rate_cummulative\n');
fprintf('    %04dx%04d    %6.2e       %6.2e       %6.2e             %6.2e\n', [2.^gridsizes; 2.^gridsizes; L2_err; L2_ratio; L2_rate; L2_rate_c]);

rate = L2_rate(end);

fprintf('\nDone testing mesh convergence of Poisson solver.\n\n');

end % function conv_test
