function plotresults(order)

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

publish_results = 0;

data_obj = importdata(['state_',int2str(order),'.txt'], ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure, trisurf(adj, nodes(:,1), nodes(:,2), state);
%hold on
%plot3(sensor(:,1),sensor(:,2),max(state)*ones(size(sensor(:,1))),'k.','MarkerSize',10)
%hold off
shading interp;
view(2);
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Computed state');
colorbar

data_obj = importdata(['control_',int2str(order),'.txt'], ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure, trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(2);
axis equal;
axis tight;
xlabel('x');
ylabel('y');
title('Computed control');
colorbar

set(findall(gcf,'-property','FontSize'),'FontSize',16); set(gcf, 'Color', 'White');% tightfig;
if publish_results
  myaa(2);
  myaa('publish');
end

end
