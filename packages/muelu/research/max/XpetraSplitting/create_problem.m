% add path to call Matrix Market routines
addpath('../MatrixMarket_routines/')

%number of nodes in x direction for each region
nregion_nodes_x = 4;

%number of nodes in y direction for each region
nregion_nodes_y = 4;

%number of intervals splitting the x dicretion of the domain
nintervals_x = 1;

%number of intervals splitting the y direction of the domain
nintervals_y = 2;

%total number of regions
nregions = nintervals_x * nintervals_y;

%Output file names
node_filename = 'node_multiregional.txt';
matrix_filename = 'A_multiregional.mm';

% routine that creates nodes structure and composite stiffness matrix
[nodes, A] = create_matrix(nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y);

%Writing node data structure and operator into a file
node_fileID = fopen(node_filename,'w');
if ( node_fileID == -1 )
 error('Cannot open file for output');
end;
fprintf(node_fileID, '%s \n', '%TotalNodes %TotalRegions');
fprintf(node_fileID, '%d \t %d\n', size(A,1), nregions);
fprintf(node_fileID, '%s\n', '%Node %Region');
fclose(node_fileID);
dlmwrite(node_filename, nodes, '-append', 'delimiter', '\t', 'precision', '%8d');

%Writing the matrix in a file with Matrix Market format
mmwriteNonOptimized(matrix_filename, A);

