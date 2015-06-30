A = laplacianfun([27, 27]);
b = (1:(27 * 27))';
matlabSoln = 0 * b; %Set 'initial guesses' to all 0
%Set up the problem using a Matlab TwoLevelFactory for Aggregates
matlabProblem = muelu('setup', A, 'xml parameter file', 'Tests/Brick/matlabParams.xml');
matlabSoln = muelu(matlabProblem, b);
