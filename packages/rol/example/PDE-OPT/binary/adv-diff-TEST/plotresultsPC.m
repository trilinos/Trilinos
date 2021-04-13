function plotresultsPC(order)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

publish_results = 0;

sensor = load('sensor.txt');
figure,
scatter(sensor(:,1),sensor(:,2),30,sensor(:,3),'filled','o')
axis equal
xlim([0,2]);
ylim([0,1])

axsize = 400;
figure('Position', [100 100 3.3*axsize 1.6*axsize])

subplot(1,2,1)
data_obj = importdata(['state_',int2str(order),'.txt'], ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
trisurf(adj, nodes(:,1), nodes(:,2), state);
hold on
plot3(sensor(:,1),sensor(:,2),max(state)*ones(size(sensor(:,1))),'k.','MarkerSize',10)
hold off
shading interp;
view(2);
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Computed state');
colorbar

n = 2^order;
subplot(1,2,2)
control = load(['control_binary_',int2str(order),'.txt']);
X       = load(['X_',int2str(order),'.txt']);
Y       = load(['Y_',int2str(order),'.txt']);
control = sparse(X+1,n-Y,control).';
imagesc(control)
axis equal;
axis tight;
box on;
xlabel('x');
ylabel('y');
title('Computed controls');
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White');% tightfig;
if publish_results
  myaa(2);
  myaa('publish');
end

end
