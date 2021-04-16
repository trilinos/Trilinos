function plotresultsC0(nx)

tag = ['_',int2str(nx),'x',int2str(nx)];
adj = load(['cell_to_node_quad',tag,'.txt']) + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load(['nodes',tag,'.txt']);  %% load node coordinates

data_obj = importdata(['control',tag,'.txt'], ' ', 2);  %% we need to skip the first two lines
ctrl = data_obj.data;

figure
xnodes = nodes(:,1);
ynodes = nodes(:,2);
patch(xnodes(adj)',ynodes(adj)',ones(size(xnodes(adj)')),'CData',ctrl,'FaceColor','flat','EdgeColor','none');
colormap(flipud(gray))
axis equal;
axis tight;
box on;
print('-depsc2',['control',tag,'.eps']);
