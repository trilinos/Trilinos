adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

publish_results = 0;

axsize = 400;
figure('Position', [100 100 3.3*axsize 1.6*axsize])

subplot(1,2,1)
data_obj = importdata('mean_state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed state mean');
colorbar

subplot(1,2,2)
data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
title('Computed controls');
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White');% tightfig;
if publish_results
  myaa(2);
  myaa('publish');
end
