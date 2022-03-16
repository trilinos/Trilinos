function vtkwrite(filename, t, p, u, c, varname)
% vtkwrite
% creates a vtk file filename.vtk containing mesh and vertex data
% supported cell types:  'Line'
%                        'Triangle'
%                        'Quadrilateral'
%                        'Tetrahedron'
%                        'Hexahedron'
%
% input: filename   vtk file destination
%                   (string)
%        t          mesh connectivity array
%                   (M x nvpc array, where M denotes the number of cells
%                    and nvpc is the number of vertices per cell )
%        p          coordinates of the mesh vertices
%                   (N x d array, where N is the number of vertices and d the dimension)
%        u          scalar values associated with each vertex
%                   (N x 1 array)
%        c          cell type
%                   (string)
%        varname    variable name
%
% example usage:
%        2d: [X,Y] = meshgrid(0:0.01:1, 0:0.01:1);
%            p = [reshape(X,prod(size(X)),1), reshape(Y,prod(size(Y)),1)];
%            t = delaunayn(p);
%            u = peaks(6*p(:,1)-3,6*p(:,2)-3);
%            vtkwrite('test2d', t, p, u, 'Triangle');
%        3d: p = rand(10, 3);
%            t = delaunayn(p);
%            u = sum(p.^2, 2);
%            vtkwrite('test3d', t, p, u, 'Tetrahedron');
%

% Extract array sizes.
[np, dim]  = size(p);  % number of mesh vertices, dimension
[nc, nvpc] = size(t);  % number of cells, number of vertices per cell

% These are VTK-specific mappings of topologies to integer designators.
ctm = containers.Map();
ctm('Line') = 3;
ctm('Triangle') = 5;
ctm('Quadrilateral') = 9;
ctm('Tetrahedron') = 10;
ctm('Hexahedron') = 12;

% Write preamble.
FID = fopen(strcat(filename, '.vtk'), 'w+');
fprintf(FID, '# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\n');
fprintf(FID, 'DATASET UNSTRUCTURED_GRID\n');

% Write points -- always pad for 3D.
fprintf(FID, 'POINTS %d float\n', np);
s = '%f %f %f \n';
P = [p zeros(np,3-dim)];
fprintf(FID, s, P');

% Write number of vertices per cell + full cell connectivity.
fprintf(FID, 'CELLS %d %d\n', nc, nc*(nvpc+1));
s = '%d ';
for k = 1:nvpc
  s = horzcat(s, {' %d'});
end
s = cell2mat(horzcat(s, {' \n'}));
fprintf(FID, s, [nvpc*ones(nc,1) t-1]');

% Write cell types.
fprintf(FID, 'CELL_TYPES %d\n', nc);
s = '%d\n';
fprintf(FID, s, ctm(c)*ones(nc,1));

% Write data.
%fprintf(FID, 'POINT_DATA %s\nSCALARS u float 1\nLOOKUP_TABLE default\n', num2str(np));
fprintf(FID, 'POINT_DATA %s\nSCALARS ', num2str(np));
fprintf(FID, varname);
fprintf(FID, ' float 1\nLOOKUP_TABLE default\n');
s = '%f\n';
fprintf(FID, s, u);

% Done.
fclose(FID);
