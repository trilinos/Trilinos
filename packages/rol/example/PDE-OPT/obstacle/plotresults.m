adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(0,90)
axis square
title('Stress potential')

data_obj = importdata('upperBound.txt', ' ', 2);  %% we need to skip the first two lines
upper = data_obj.data;
figure
trisurf(adj, nodes(:,1), nodes(:,2), upper);
shading interp;
view(0,90)
axis square
title('Upper bounds')

figure
trisurf(adj, nodes(:,1), nodes(:,2), (state>=upper).*ones(size(state)));
shading interp;
view(0,90)
colormap bone
colormap(flipud(colormap))
axis square
title('Active set')