adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

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
quiver(nodes(:,1), nodes(:,2), Ux, Uy);
axis equal
axis tight

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
data_obj = importdata('map_control.txt', ' ', 9);  %% we need to skip the first 9 lines
map_control = data_obj.data;
map_control = map_control(1:2:end)+1;
[tmp, control_permute] = sort(map_control);
control = control(control_permute);  %% we need to permute the control according to parallel maps

% Extract x-velocity
Zx = control(1:3:N);
% Extract y-velocity
Zy = control(2:3:N);
% Extract pressure
P  = control(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zx);
shading interp;
view(0,90)
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zy);
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
quiver(nodes(:,1), nodes(:,2), Zx, Zy);
axis equal
axis tight

figure,
idx = find((nodes(:,1) == 1) .* (nodes(:,2) <= 0.5));
ZX = zeros(size(Zx)); ZX(idx) = Zx(idx);
ZY = zeros(size(Zy)); ZY(idx) = Zy(idx);
quiver(nodes(:,1), nodes(:,2), ZX, ZY);
axis equal
axis tight
xlim([0.6,1.2]);
ylim([-0.2,0.7]);

figure,
idx = find((nodes(:,1) == 1) .* (nodes(:,2) <= 0.5));
UX = zeros(size(Ux)); UX(idx) = Ux(idx);
UY = zeros(size(Uy)); UY(idx) = Uy(idx);
quiver(nodes(:,1), nodes(:,2), UX, UY);
axis equal
axis tight
xlim([0.6,1.2]);
ylim([-0.2,0.7]);
