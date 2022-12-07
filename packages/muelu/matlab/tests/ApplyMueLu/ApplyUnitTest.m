addpath '../../bin' '../Common'

try
  [A, coords] = laplacianfun([90, 90]);
  b = (1:(8100))';
  % Set up the problem using a Matlab TwoLevelFactory for Aggregates
  mueluProblem = muelu('setup', A, coords, 'xml parameter file', 'mueluParams.xml');
  % Create initial guess (of zeros)
  x = zeros(8100, 1);
  % Compute initial residual norm
  initErr = norm(A*x - b);
  % Apply MueLu preconditioner for 10 iterations
  for i = 0:25
    r = A*x - b;
    fprintf('Residual norm at iter %d: %f\n', i, norm(r) / initErr);
    update = muelu('apply', mueluProblem, r);
    x = x - update;
  end
  muelu('cleanup', mueluProblem);
  % Check that final residual norm is small
  if norm(A*x - b) / initErr < 1e-2
    disp('Test passed.');
    exit(0);
  else
    disp('Apply test failed: residual decreased less than expected.');
    exit(-1);
  end
catch me
  disp('Test failed, on exception:');
  disp(getReport(me));
  exit(-2)
end
