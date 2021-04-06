adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

figure('Renderer','painters','Position',[10 10 1800 600])

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
subplot(1,3,1)
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(0,90)
axis square
title('Temperature','FontSize',20)
set(gca,'FontSize',16)

data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
if (length(control) < length(state))
  subplot(1,3,2)
  xnodes = nodes(:,1);
  ynodes = nodes(:,2);
  patch(xnodes(adj)',ynodes(adj)',ones(size(xnodes(adj)')),'CData',control,'FaceColor','flat','EdgeColor','none');
  axis equal;
  axis tight;
  box on;
  axis square
else
  subplot(1,3,2)
  trisurf(adj, nodes(:,1), nodes(:,2), control);
  shading interp;
  view(0,90)
  axis square
end
title('Density','FontSize',20)
set(gca,'FontSize',16)

data_obj = importdata('filtered_density.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
subplot(1,3,3)
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
view(0,90)
axis square
title('Filtered Density','FontSize',20)
set(gca,'FontSize',16)
