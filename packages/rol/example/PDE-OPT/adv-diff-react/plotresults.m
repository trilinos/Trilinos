adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

figure()
data_obj = importdata('weights.txt', ' ', 2);  %% we need to skip the first two lines
weights = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), weights);
shading flat;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Sensors');
colormap summer
set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White'); tightfig;
myaa(2); myaa('publish');

axsize = 400;
figure('Position', [100 100 3*axsize 1.6*axsize])

subplot(2,3,4)
data_obj = importdata('data.txt', ' ', 2);  %% we need to skip the first two lines
data = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), data);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Data (true state)');

subplot(2,3,1)
data_obj = importdata('sources.txt', ' ', 2);  %% we need to skip the first two lines
sources = data_obj.data;
% normalize source data to 1
maxsource = max(abs(sources));
sources = sources * (1/maxsource);
trisurf(adj, nodes(:,1), nodes(:,2), sources);
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

subplot(2,3,2)
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
trisurf(adj, nodes(:,1), nodes(:,2), sources-control);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Source error (rel)');
colorbar

subplot(2,3,6)
trisurf(adj, nodes(:,1), nodes(:,2), (state-data)/max(abs(data)));
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('State error (rel)');
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White'); tightfig;
myaa(2); myaa('publish');
