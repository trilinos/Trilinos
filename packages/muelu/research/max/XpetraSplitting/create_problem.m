addpath('../MatrixMarket_routines/')

nregion_nodes_x = 4;
nregion_nodes_y = 4;
nintervals_x = 2;
nintervals_y = 2;


[nodes, A] = create_matrix(nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y);


dlmwrite('nodes_multiregional.txt', nodes, '\t');
mmwriteNonOptimized('A_multiregional.mtx', A);

