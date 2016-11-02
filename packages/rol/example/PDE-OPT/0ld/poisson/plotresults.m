adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
subplot(2,2,1)
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(0,90)
axis square
title('State (Reduced-space)')

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
subplot(2,2,2)
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(0,90)
axis square
title('Control (Reduced-space)')

data_obj = importdata('stateFS.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
subplot(2,2,3)
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(0,90)
axis square
title('State (Full-space)')

data_obj = importdata('controlFS.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
subplot(2,2,4)
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(0,90)
axis square
title('Control (Full-space)')
