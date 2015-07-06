A = laplacianfun([90, 90]);
b = (1:(8100))';
%Set up the problem using a Matlab TwoLevelFactory for Aggregates
matlabProblem = muelu('setup', A, 'xml parameter file', 'Tests/Brick/matlabParams.xml');
matlabP = muelu('get', matlabProblem, 1, 'P');
disp('MATLAB brick aggregation test passed by running to completion.');
