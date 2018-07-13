adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
stateX = data_obj.data(1:2:end);
stateY = data_obj.data(2:2:end);
subplot(1,2,1)
trisurf(adj, nodes(:,1), nodes(:,2), stateX);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed displacement X');
subplot(1,2,2)
trisurf(adj, nodes(:,1), nodes(:,2), stateY);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed displacement Y');
