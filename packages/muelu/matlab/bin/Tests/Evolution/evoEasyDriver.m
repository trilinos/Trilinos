%Driver script, works as a unit test
%Will force MATLAB to return with exit code depending on comparison
%of tentative prolongator from matlab function and muelu
addpath('../..')

A = laplacianfun([80, 80]);
matlabProblem = muelu('setup', A, 'xml parameter file', 'matlabParamsEasy.xml');
