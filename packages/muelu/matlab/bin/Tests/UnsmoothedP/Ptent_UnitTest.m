%Driver script, works as a unit test
%Will force MATLAB to return with exit code depending on comparison
%of tentative prolongator from matlab function and muelu

%Create a 2D laplacian problem (80x80 gets 2 MueLu levels)
A = laplacianfun([60, 60]);
%Note: In both of these problems, set reuse type to full so that P is not discarded
%Set up the problem normally in MueLu through MueMex
mueluProblem = muelu('setup', A, 'reuse: type', 'full');
%Set up the problem using a Matlab TwoLevelFactory for Ptent
matlabProblem = muelu('setup', A, 'xml parameter file', 'Tests/UnsmoothedP/params.xml');
mueluP = muelu('get', mueluProblem, 1, 'P');
matlabP = muelu('get', matlabProblem, 1, 'P');
%Now compare the matrices to near machine precision for equality
areEqual = true(1)
for index = 1:numel(mueluP)
  %Use machine precision minus a few digits
  precision = eps(mueluP(index)) * 1000;
  if(abs(mueluP(index) - matlabP(index)) > precision)
    areEqual = false(1)
    break
  end
end
if areEqual
  disp('Test passed, MueLu tentative prolongator is correct.');
else
  disp('Test failed, MueLu tentative prolongator did not match MATLAB one.');
end
