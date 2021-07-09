%Driver script for reentrant mumex test
%Sets up a problem where for level 1, P is created by making a fresh hierarchy
%from within MATLAB and modifying it.

addpath('../../bin')

try
  A = laplacianfun([120, 120]);
  matlabProblem = muelu('setup', A, 'xml parameter file', 'matlabParams.xml');
  muelu('cleanup');
  disp('Test passed by running to completion with reentrant call.');
  exit(0);
catch me
  fprintf('Test failed with exception "%s"\n', me.message);
  exit(-2);
end
