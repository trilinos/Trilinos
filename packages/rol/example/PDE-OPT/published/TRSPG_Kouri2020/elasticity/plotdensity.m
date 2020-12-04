function plotdensity

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('density.txt', ' ', 2);  %% we need to skip the first two lines
dens = data_obj.data(1:2:end,:);

figure
xnodes = nodes(:,1);
ynodes = nodes(:,2);
celldens = dens(adj);
avgcelldens = sum(celldens,2)/4;
%maxcelldens = max(celldens,[],2);
patch(xnodes(adj)',ynodes(adj)',ones(size(xnodes(adj)')),'CData',avgcelldens,'FaceColor','flat','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
