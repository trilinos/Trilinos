addpath '../../bin' '../Common'

try
  [A, coords] = laplacianfun([90, 90]);
  b = (1:(8100))';
  % Set up the problem using a Matlab TwoLevelFactory for Aggregates
  matlabProblem = muelu('setup', A, coords, 'xml parameter file', 'matlabParams.xml');
  mueluProblem = muelu('setup', A, coords, 'xml parameter file', 'mueluParams.xml');
  matlabP = muelu('get', matlabProblem, 1, 'P');
  mueluP = muelu('get', mueluProblem, 1, 'P');
  muelu('cleanup');
  diff = nonzeros(matlabP - mueluP);
  passed = true(1);
  for i = 1:numel(diff)
    elem = abs(diff(i));
    if elem > 1e-12
      passed = false(1);
      break
    end
  end
  if passed
    disp('Brick aggregation test passed.');
    exit(0);
  else
    disp('Brick aggregation test failed.');
    exit(-1);
  end

catch me
  disp('Test failed, on exception:');
  disp(getReport(me));
  exit(-2)
end
