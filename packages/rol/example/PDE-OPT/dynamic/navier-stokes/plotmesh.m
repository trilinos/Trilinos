adj   = load('cell_to_node_quad.txt') + 1; %% load node adjacency, increment by 1 for 1-based indexing
nodes = load('nodes.txt');                 %% load node coordinates

x = nodes(:,1);
y = nodes(:,2);
d = adj(:,[1 2 3 4 1]).';

figure,
plot(x(d), y(d), 'k')
axis('equal','tight')
print('mesh.eps','-depsc2','-r800');
