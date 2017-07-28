% add path to call Matrix Market routines
addpath('../MatrixMarket_routines/')

%number of nodes in x direction for each region
nregion_nodes_x = 4;

%number of nodes in y direction for each region
nregion_nodes_y = 4;

%number of intervals splitting the x dicretion of the domain
nintervals_x = 2;

%number of intervals splitting the y direction of the domain
nintervals_y = 2;

%Output file names
node_filename = 'nodes_multiregional.txt';
matrix_filename = 'A_multiregional.txt';

% routine that creates nodes structure and composite stiffness matrix
[nodes, A] = create_matrix(nregion_nodes_x, nregion_nodes_y, nintervals_x, nintervals_y);

%Writing node data structure and operator into a file
file_opened = fopen([node_filename],'w');
if ( file_opened == -1 )
 error('Cannot open file for output');
end;
fprintf(node_filename, '%TotalNodes \t %TotalRegions\n');
fprintf(node_filename, '%d \t %d\n', size(A,1), nregion_nodes_x*nregion_nodes_y);
fprintf(node_filename, '%Node \t %lRegion\n');
fclose(node_filename);
dlmwrite(filename, nodes, '-append', '\t');

%Writing the matrix in a file with Matrix Market format
mmwriteNonOptimized(matrix_filename, A);

