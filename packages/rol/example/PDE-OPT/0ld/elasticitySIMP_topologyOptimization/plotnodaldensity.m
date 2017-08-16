adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('density_nodal.txt', ' ', 2);  %% we need to skip the first two lines
dens = data_obj.data;
figure
trisurf(adj, nodes(:,1), nodes(:,2), dens);
shading interp;
view(2);
axis tight;
xlabel('x');
ylabel('y');
title('Nodal density');
