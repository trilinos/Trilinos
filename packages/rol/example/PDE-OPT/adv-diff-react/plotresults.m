adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

figure('Position', [100 100 1000 400])

subplot(1,2,1)
data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');

subplot(1,2,2)
data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
