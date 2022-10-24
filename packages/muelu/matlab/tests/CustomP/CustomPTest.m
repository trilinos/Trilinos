addpath('../../bin')

try
  A = laplacianfun([50, 50]);
  KeepNodes = int32([5, 9, 23, 63, 120, 121, 150, 200, 300, 634, 895]); %These nodes will become singleton aggregates (in P)
  %Set up the problem using a Matlab TwoLevelFactory for P
  matlabProblem = muelu('setup', A, 'coarse: max size', 25, 'xml parameter file', 'matlabParams.xml', 'level 0', {'OrdinalVector KeepNodes', KeepNodes});
  muelu('cleanup');
  disp('Custom aggregation test passed by running to completion.');
  exit(0);
catch me
  disp('Test failed, on exception:');
  disp(getReport(me));
  exit(-2);
end
