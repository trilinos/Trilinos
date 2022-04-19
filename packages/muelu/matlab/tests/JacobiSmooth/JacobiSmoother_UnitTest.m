%Driver script, works as a unit test
%Will force MATLAB to return with exit code depending on comparison
%of solution vectors from matlab function and muelu
addpath('../../bin')

try
  %Create a small problem with only one level
  A = laplacianfun([80, 80]);
  b = (1:6400)';
  matlabSoln = 0 * b; %Set 'initial guesses' to all 0
  mueluSoln = matlabSoln;
  %Note: In both of these problems, set reuse type to full so that P is not discarded
  mueluProblem = muelu('setup', A, 'xml parameter file', 'mueluParams.xml');
  matlabProblem = muelu('setup', A, 'xml parameter file', 'matlabParams.xml');
  mueluSoln = muelu(mueluProblem, b);
  matlabSoln = muelu(matlabProblem, b);
  %Compare the solutions to near machine precision for equality
  areEqual = true(1);
  for index = 1:numel(matlabSoln)
    %Use machine precision minus a few digits
    precision = eps(mueluSoln(index)) * 1000;
    if(abs(mueluSoln(index) - matlabSoln(index)) > precision)
      areEqual = false(1);
      break
    end
  end
  muelu('cleanup');
  if areEqual
    disp('Test passed, MueLu produced same solution as gold standard in MATLAB.');
    exit(0);
  else
    disp('Test failed, MueLu''s solution did not match MATLAB''s.');
    exit(-1);
  end

catch me
  disp('Test failed, on exception:');
  disp(getReport(me));
  exit(-2)
end
