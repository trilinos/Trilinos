function plotnodaldensity(nx,ny)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
dens = data_obj.data(1:2:end,:);
figure
trisurf(adj, nodes(:,1), nodes(:,2), dens);
shading flat;
view(2);
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Nodal density');

figure
v1 = 0.25;
v2 = 0.25;
V  = [v1,v2];
contourf(reshape(nodes(:,1),nx+1,ny+1), reshape(nodes(:,2),nx+1,ny+1), reshape(dens,nx+1,ny+1), V);
colormap([0 0 0]);
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Nodal density');

figure
quadsurf(length(nodes), nodes.', nx*ny, adj.', dens)
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Nodal density');