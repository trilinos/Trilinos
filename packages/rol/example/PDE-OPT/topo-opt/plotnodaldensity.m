function plotnodaldensity(nx,ny)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
dens = data_obj.data(1:2:end,:);
figure
patch('Faces',adj,'Vertices',nodes,'FaceVertexCData',dens,'FaceColor','interp','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
xlabel('x');,
ylabel('y');
title('Interpolated nodal density');

figure
xnodes = nodes(:,1);
ynodes = nodes(:,2);
celldens = dens(adj);
avgcelldens = sum(celldens,2)/4;
%maxcelldens = max(celldens,[],2);
patch(xnodes(adj)',ynodes(adj)',ones(size(xnodes(adj)')),'CData',avgcelldens,'FaceColor','flat','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Averaged nodal density');

figure
v1 = 0.5;
v2 = 0.5;
V  = [v1,v2];
contourf(reshape(nodes(:,1),nx+1,ny+1), reshape(nodes(:,2),nx+1,ny+1), reshape(dens,nx+1,ny+1), V);
colormap([0 0 0]);
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Contoured nodal density');