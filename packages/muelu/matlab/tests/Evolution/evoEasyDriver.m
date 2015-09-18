%Driver script, works as a unit test
%Will force MATLAB to return with exit code depending on comparison
%of tentative prolongator from matlab function and muelu
addpath('../../bin')

%filename = 'Aiso.mtx';
%filename = 'Aani.mtx';
filename = 'Arot.mtx';
fprintf('Reading %s\n',filename);
A = mmread(filename);
matlabProblem = muelu('setup', A, 'xml parameter file', 'matlabParamsEasy.xml');
