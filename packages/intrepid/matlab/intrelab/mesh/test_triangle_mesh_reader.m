%
% AUTHOR: Denis Ridzal
%         Computational and Applied Mathematics
%         Rice University
%

% get point, triangle and edge information from file 
% '../mesh2_datagrid_airport/airport.1'
filename = '../mesh2_data/grid_airport/airport.1';

fprintf(1,' read mesh data from  %s \n', filename)
[mesh] = triangle_mesh_reader( filename );

figure(1)
plot_mesh2(mesh);
axis('tight','equal')

% refine the mesh
% NOTE: Currently, this currently only works properly for 
%       triangles wiht 3 nodes. If triangles have 6 nodes,
%       only mesh.t(:,1:3) is used

fprintf(1,' refine the mesh  \n')
[mesh1,I1]=reg_mesh_refine2(mesh);
% plot the mesh
figure(2)
plot_mesh2(mesh1);
axis('tight','equal')

fprintf(1,' refine the mesh again \n')
[mesh2,I2]=reg_mesh_refine2(mesh1);
% plot the mesh
figure(3)
plot_mesh2(mesh2);
axis('tight','equal')
