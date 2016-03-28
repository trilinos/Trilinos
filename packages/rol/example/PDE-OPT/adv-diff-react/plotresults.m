adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

axsize = 400;
figure('Position', [100 100 3*axsize 1.6*axsize])

subplot(2,3,2)
data_obj = importdata('data.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Data (true state)');

subplot(2,3,1)
data_obj = importdata('sources.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('True sources');

subplot(2,3,5)
data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed state');

subplot(2,3,4)
data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed sources');

subplot(2,3,3)
data_obj = importdata('weights.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading flat;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Sensors');
