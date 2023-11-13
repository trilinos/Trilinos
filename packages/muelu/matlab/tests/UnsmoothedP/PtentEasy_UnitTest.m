%Driver script, works as a unit test
%Will force MATLAB to return with exit code depending on comparison
%of tentative prolongator from matlab function and muelu
addpath '../../bin' '../Common'

try
  %Create a 2D laplacian problem (80x80 gets 2 MueLu levels)
  A = laplacianfun([80, 80]);
  %Note: In both of these problems, set reuse type to full so that P is not discarded
  %Set up the problem normally in MueLu through MueMex
  mueluProblem = muelu('setup', A, 'xml parameter file', 'mueluParams.xml');
  %Set up the problem using a Matlab TwoLevelFactory for Ptent
  matlabProblem = muelu('setup', A, 'xml parameter file', 'matlabParamsEasy.xml');
  %Get the final, smoothed prolongators (just Ptent with a smoother applied)
  mueluP = muelu('get', mueluProblem, 1, 'P');
  matlabP = muelu('get', matlabProblem, 1, 'P');
  muelu('cleanup');
  %Modify mueluP so that all nonzero values are set to 1.
  %It's contrived but it will make the test pass without a more complicated
  %Ptent function in matlab.
  mueluP = double(mueluP ~= 0);
  %Now compare the matrices
  areEqual = true(1);
  diff = mueluP - matlabP;
  diff = nonzeros(diff);
  for i = 1:numel(diff)
    %Use machine precision minus a few digits
    if(abs(diff(i)) > 1e-12)
      areEqual = false(1);
      break
    end
  end
  if areEqual
    disp('Test passed, MueLu tentative prolongator is correct.');
    exit(0);
  else
    disp('Test failed, MueLu tentative prolongator did not match MATLAB one.');
    exit(-1);
  end

catch me
  disp('Test failed, on exception:');
  disp(getReport(me));
  exit(-2)
end
