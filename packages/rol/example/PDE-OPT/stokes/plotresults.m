adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

N = 3*length(nodes);
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
title('X-Velocity')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
title('Y-Velocity')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(Ux.^2 + Uy.^2));
shading interp;
view(0,90)
title('Velocity Magnitude')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
view(0,90)
title('Pressure')
axis equal
axis tight

figure,
quiver(nodes(:,1), nodes(:,2), Ux, Uy);
title('Velocity')
axis equal
axis tight
