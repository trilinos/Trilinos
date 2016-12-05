adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;

N = 4*length(nodes);
% Extract x-velocity
Ux = state(1:4:N);
% Extract y-velocity
Uy = state(2:4:N);
% Extract pressure
P  = state(3:4:N);
% Extract temperature
T  = state(4:4:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), sqrt(Ux.^2 + Uy.^2));
shading interp;
view(0,90)
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
view(0,90)
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), T);
shading interp;
view(0,90)
axis equal
axis tight

figure,
quiver(nodes(:,1), nodes(:,2), Ux, Uy);
axis equal
axis tight

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

% Extract x-velocity
Zx = control(1:4:N);
% Extract y-velocity
Zy = control(2:4:N);
% Extract pressure
P  = control(3:4:N);
% Extract Temperature
T  = control(4:4:N);

Xmin = min(nodes(:,1));
Xmax = max(nodes(:,1));

figure,
idx = find(nodes(:,2) == 0.0);
plot(nodes(idx,1), T(idx));
xlim([Xmin,Xmax]);

figure,
idx = find(nodes(:,2) == 1.0);
plot(nodes(idx,1), T(idx));
xlim([Xmin,Xmax]);
