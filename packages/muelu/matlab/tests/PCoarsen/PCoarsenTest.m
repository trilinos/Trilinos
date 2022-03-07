addpath('../../bin')

try
  A = laplacianfun([200],-1);
  
  % Setup an element-to-node map for a 1D Laplacian (yay)
  n=size(A,1);
  e2n=int32([[n,1:n-1];[1:n]]');
  matlabProblem = muelu('setup', A, 'coarse: max size', 25, ...
                        'verbosity','high',...
                        'multigrid algorithm','pcoarsen',...
                        'pcoarsen: schedule','{1,1}',...
                        'pcoarsen: element','hgrad_line_c',...
                        'level 0',{'pcoarsen: element to node map',e2n});
  muelu('cleanup');
  disp('P-coarsening test passed by running to completion.');
  exit(0);
catch me
  disp('Test failed, on exception:');
  disp(getReport(me));
  exit(-2);
end
